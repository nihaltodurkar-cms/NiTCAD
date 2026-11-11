// The TVB-31 parameter register (Unit 55, TVR-3; DECISIONS 15.1 (b), 15.9, 15.10). Data only, base-free.
//
// One row per model parameter NiTCAD ships. `source_value` is typed from the source as printed (digits and exponent
// only; units in `unit`), never from NiTCAD's code or structs (§2.2). Where no source has been read, it is empty and
// the grade is L: the table check of that row cannot pass and is reported as "provenance insufficient" (MI-Q5).
// `binding` names the NiTCAD field the benchmark compares, as text; the reference never reads that field.
//
// Kinds:
//   printed  - the source prints this value; the shipped value must equal its correctly rounded parse (15.9);
//   law_300K - the source prints a temperature law; `source_value` is the law's value at 300 K computed here (PD,
//              shown in `note`), and the shipped value is judged to the source's printed digits of the law's
//              coefficients (half a unit in the last place, the NPARC rule of 11.5);
//   choice   - tag C: a NiTCAD or owner choice, pinned to the ruling's value, not verified.
//
// Grades (§3.3 and MI-Q5): PA primary read; PB two or more independent reproductions agree; PC one reproduction;
// PD derived here; L legacy value, no source read; C choice.
#pragma once

#include <array>
#include <string_view>

namespace NiTCAD::verify::reference::model_parameters {

enum class Grade { PA, PB, PC, PD, L, C };
enum class Kind { printed, law_300K, choice };

struct Row {
    std::string_view entry;         // TVB-31 entry, "31.06"
    std::string_view binding;       // NiTCAD field compared, e.g. "legacy_silicon.electron_mobility.mu_min"
    std::string_view source_value;  // as printed (or the law's 300 K value for law_300K); empty for grade L
    std::string_view unit;
    Kind kind;
    Grade grade;
    std::string_view source;        // the source and its locator (equation, table)
    std::string_view note;
};

// Abbreviations of the reproductions read in step 1 (DECISIONS 15.10):
//   Palankovski = V. Palankovski, PhD thesis, TU Wien 2000; Allpix = Allpix Squared 3.3 manual; Genius = Genius
//   TCAD Open, src/material/Si; Xyce = Xyce reference guide tables and N_DEV_Material*.C; Tidy3D = Tidy3D docs.
inline constexpr std::array rows{
    // 31.01 Varshni, legacy silicon.
    Row{"31.01", "legacy_silicon.Eg0_eV", "1.17", "eV", Kind::printed, Grade::PB,
        "Thurmond 1975 as reproduced by Xyce bandgapTbl.tex; Palankovski Table 3.8 (reported value)",
        "Palankovski's fitted value is 1.1695"},
    Row{"31.01", "legacy_silicon.varshni_alpha_eV_per_K", "4.73e-4", "eV/K", Kind::printed, Grade::PB,
        "Xyce bandgapTbl.tex; Palankovski Table 3.8", ""},
    Row{"31.01", "legacy_silicon.varshni_beta_K", "636", "K", Kind::printed, Grade::PB,
        "Xyce bandgapTbl.tex; Palankovski Table 3.8", ""},
    Row{"31.01", "legacy_silicon.Nc300", "2.86e19", "cm^-3", Kind::printed, Grade::PB,
        "Green, J. Appl. Phys. 67, 2944 (1990), as Genius Si_band_default.cc NC300; Palankovski Table 3.20 sum",
        ""},
    Row{"31.01", "legacy_silicon.Nv300", "3.10e19", "cm^-3", Kind::printed, Grade::PB,
        "Green 1990, as Genius NV300; Palankovski Table 3.20 sum", "Xyce ships 2.66e19 (\"from Larry\")"},
    Row{"31.01", "legacy_silicon.eps_r", "", "", Kind::printed, Grade::L, "not read", ""},
    Row{"31.01", "legacy_silicon.electron_affinity_eV", "", "eV", Kind::printed, Grade::L, "not read", ""},
    // 31.02, 31.03 Couderc et al. (2014), read in Unit 35 (DECISIONS Unit 35 table).
    Row{"31.02", "silicon.Eg0_eV", "1.17", "eV", Kind::printed, Grade::PA,
        "Couderc, Amara, Lemiti, J. Appl. Phys. 115, 093705 (2014), Table I", ""},
    Row{"31.02", "silicon.passler_gap.alpha_eV_per_K", "3.23e-4", "eV/K", Kind::printed, Grade::PA,
        "Couderc 2014, Table I", ""},
    Row{"31.02", "silicon.passler_gap.Theta_K", "446", "K", Kind::printed, Grade::PA, "Couderc 2014, Table I", ""},
    Row{"31.02", "silicon.passler_gap.Delta", "0.51", "", Kind::printed, Grade::PA, "Couderc 2014, Table I", ""},
    Row{"31.03", "silicon.dos_mass.conduction[0]", "-4.609e-10", "K^-3", Kind::printed, Grade::PA,
        "Couderc 2014, eq (18)", "signs fixed by the paper's cross-checks (Unit 35)"},
    Row{"31.03", "silicon.dos_mass.conduction[1]", "6.753e-7", "K^-2", Kind::printed, Grade::PA,
        "Couderc 2014, eq (18)", ""},
    Row{"31.03", "silicon.dos_mass.conduction[2]", "-1.312e-5", "K^-1", Kind::printed, Grade::PA,
        "Couderc 2014, eq (18)", ""},
    Row{"31.03", "silicon.dos_mass.conduction[3]", "1.094", "", Kind::printed, Grade::PA, "Couderc 2014, eq (18)",
        ""},
    Row{"31.03", "silicon.dos_mass.valence[0]", "2.525e-9", "K^-3", Kind::printed, Grade::PA,
        "Couderc 2014, eq (19)", ""},
    Row{"31.03", "silicon.dos_mass.valence[1]", "-4.689e-6", "K^-2", Kind::printed, Grade::PA,
        "Couderc 2014, eq (19)", ""},
    Row{"31.03", "silicon.dos_mass.valence[2]", "3.376e-3", "K^-1", Kind::printed, Grade::PA,
        "Couderc 2014, eq (19)", ""},
    Row{"31.03", "silicon.dos_mass.valence[3]", "0.3426", "", Kind::printed, Grade::PA, "Couderc 2014, eq (19)",
        ""},
    // 31.05 band-gap narrowing (F1).
    Row{"31.05", "legacy_silicon.bandgap_narrowing.E0_eV", "6.92e-3", "eV", Kind::printed, Grade::PB,
        "Klaassen, Slotboom, de Graaff, Solid-State Electron. 35, 125 (1992), as Tidy3D and Genius V0_BGN",
        "F1: the 1976 Slotboom-de Graaff value is 9e-3 (Palankovski Table 3.15); the header cites Slotboom"},
    Row{"31.05", "legacy_silicon.bandgap_narrowing.N0", "1.3e17", "cm^-3", Kind::printed, Grade::PB,
        "Klaassen et al. 1992, as Tidy3D and Genius N0_BGN", "F1: 1976 value 1e17"},
    // 31.06 Caughey-Thomas (F3): no source read prints these sets.
    Row{"31.06", "legacy_silicon.electron_mobility.mu_min", "", "cm^2/(V s)", Kind::printed, Grade::L,
        "not found (Palankovski, Xyce, Genius, Allpix read)", "F3"},
    Row{"31.06", "legacy_silicon.electron_mobility.mu_max", "", "cm^2/(V s)", Kind::printed, Grade::L, "not found",
        "F3"},
    Row{"31.06", "legacy_silicon.electron_mobility.N_ref", "", "cm^-3", Kind::printed, Grade::L, "not found", "F3"},
    Row{"31.06", "legacy_silicon.electron_mobility.alpha", "", "", Kind::printed, Grade::L, "not found", "F3"},
    Row{"31.06", "legacy_silicon.electron_mobility.T_exponent", "-2.33", "", Kind::printed, Grade::PB,
        "Arora, Hauser, Roulston, IEEE TED 29, 292 (1982), Table 1, as Allpix and Xyce aroraTbl.tex (exn2)",
        "F3: Arora applies it to mu_0 in its own form, not to a Caughey-Thomas mu_max"},
    Row{"31.06", "legacy_silicon.hole_mobility.mu_min", "", "cm^2/(V s)", Kind::printed, Grade::L, "not found", "F3"},
    Row{"31.06", "legacy_silicon.hole_mobility.mu_max", "", "cm^2/(V s)", Kind::printed, Grade::L, "not found", "F3"},
    Row{"31.06", "legacy_silicon.hole_mobility.N_ref", "", "cm^-3", Kind::printed, Grade::L, "not found", "F3"},
    Row{"31.06", "legacy_silicon.hole_mobility.alpha", "", "", Kind::printed, Grade::L, "not found", "F3"},
    Row{"31.06", "legacy_silicon.hole_mobility.T_exponent", "-2.23", "", Kind::printed, Grade::PB,
        "Arora et al. 1982, Table 1, as Allpix and Xyce aroraTbl.tex (exp2)", "F3, as the electron row"},
    // 31.07 Masetti (Unit 34/38 reads, and Allpix in step 1).
    Row{"31.07", "masetti_boron_silicon.mu_min1", "44.9", "cm^2/(V s)", Kind::printed, Grade::PB,
        "Masetti, Severi, Solmi, IEEE TED 30, 764 (1983), Table I, as Allpix and Unit 34's reproductions", ""},
    Row{"31.07", "masetti_boron_silicon.mu_min2", "0", "cm^2/(V s)", Kind::printed, Grade::PB,
        "Masetti 1983 eq (4) as Allpix prints it (no mu_min2 term in the hole form)", ""},
    Row{"31.07", "masetti_boron_silicon.mu_max", "470.5", "cm^2/(V s)", Kind::printed, Grade::PB, "as above", ""},
    Row{"31.07", "masetti_boron_silicon.mu_1", "29.0", "cm^2/(V s)", Kind::printed, Grade::PB, "as above", ""},
    Row{"31.07", "masetti_boron_silicon.P_c", "9.23e16", "cm^-3", Kind::printed, Grade::PB, "as above", ""},
    Row{"31.07", "masetti_boron_silicon.C_r", "2.23e17", "cm^-3", Kind::printed, Grade::PB, "as above", ""},
    Row{"31.07", "masetti_boron_silicon.C_s", "6.1e20", "cm^-3", Kind::printed, Grade::PB, "as above", ""},
    Row{"31.07", "masetti_boron_silicon.alpha", "0.719", "", Kind::printed, Grade::PB, "as above", ""},
    Row{"31.07", "masetti_boron_silicon.beta", "2.0", "", Kind::printed, Grade::PB, "as above", ""},
    Row{"31.07", "masetti_boron_silicon.T_exponent", "-2.2", "", Kind::choice, Grade::C,
        "the reproductions' Hamburg-law scaling, not Masetti's (F6)", ""},
    Row{"31.07", "masetti_phosphorus_silicon.mu_min1", "68.5", "cm^2/(V s)", Kind::printed, Grade::PB,
        "Masetti 1983 Table I, as Allpix and Palankovski Table 3.24", ""},
    Row{"31.07", "masetti_phosphorus_silicon.mu_min2", "68.5", "cm^2/(V s)", Kind::printed, Grade::PB,
        "Masetti 1983 eq (1), Table I, as Allpix (mu_0,e) and Palankovski Table 3.24", ""},
    Row{"31.07", "masetti_phosphorus_silicon.mu_max", "1414", "cm^2/(V s)", Kind::printed, Grade::PB,
        "Masetti 1983 Table I, as Allpix and Unit 38's reproduction", "raised from PC by Allpix (15.10, F6)"},
    Row{"31.07", "masetti_phosphorus_silicon.mu_1", "56.1", "cm^2/(V s)", Kind::printed, Grade::PB, "as above", ""},
    Row{"31.07", "masetti_phosphorus_silicon.C_r", "9.20e16", "cm^-3", Kind::printed, Grade::PB, "as above", ""},
    Row{"31.07", "masetti_phosphorus_silicon.C_s", "3.41e20", "cm^-3", Kind::printed, Grade::PB,
        "Allpix; Palankovski prints 3.4e20", ""},
    Row{"31.07", "masetti_phosphorus_silicon.alpha", "0.711", "", Kind::printed, Grade::PB, "as above", ""},
    Row{"31.07", "masetti_phosphorus_silicon.beta", "1.98", "", Kind::printed, Grade::PB,
        "Allpix; Palankovski prints -1.98", ""},
    Row{"31.07", "masetti_phosphorus_silicon.T_exponent", "-2.5", "", Kind::choice, Grade::C, "as boron (F6)", ""},
    // 31.08 Canali (F2).
    Row{"31.08", "legacy_silicon.electron_saturation.v_sat_cm_s", "1.0705e7", "cm/s", Kind::law_300K, Grade::PB,
        "Jacoboni et al., Solid-State Electron. 20, 77 (1977), Table 5, v_m = 1.53e9 T^-0.87, as Allpix",
        "1.53e9 * 300^-0.87 = 1.07052e7 (step 1 had 1.0704e7; corrected in step 4, 15.14); Canali's own eq 2a, 1.43e9 T^-0.87, gives 1.0005e7 (F2)"},
    Row{"31.08", "legacy_silicon.hole_saturation.v_sat_cm_s", "8.345e6", "cm/s", Kind::law_300K, Grade::PB,
        "Canali 1975 / Jacoboni 1977 Table 5, v_m = 1.62e8 T^-0.52, as Allpix",
        "1.62e8 * 300^-0.52 = 8.345e6; shipped 8.37e6 (F2)"},
    Row{"31.08", "legacy_silicon.electron_saturation.beta", "1.109", "", Kind::law_300K, Grade::PB,
        "Canali/Jacoboni beta = 2.57e-2 T^0.66, as Allpix", "shipped 2 (F2; Palankovski Table 3.27 and Genius ship "
        "2 uncited)"},
    Row{"31.08", "legacy_silicon.hole_saturation.beta", "1.213", "", Kind::law_300K, Grade::PB,
        "Canali/Jacoboni beta = 0.46 T^0.17, as Allpix", "shipped 1 (F2)"},
    // 31.09 Lombardi (Unit 46: Xyce, Genius; COMSOL as extracted).
    Row{"31.09", "lombardi_silicon_electrons.B", "4.75e7", "cm/s", Kind::printed, Grade::PB,
        "Lombardi et al., IEEE TCAD 7, 1164 (1988), as Xyce lombardiTbl.tex and Genius (Unit 46)", ""},
    Row{"31.09", "lombardi_silicon_electrons.C", "1.74e5", "", Kind::printed, Grade::PB, "as above", ""},
    Row{"31.09", "lombardi_silicon_electrons.lambda", "0.125", "", Kind::printed, Grade::PB, "as above", ""},
    Row{"31.09", "lombardi_silicon_electrons.delta", "5.82e14", "V/s", Kind::printed, Grade::PB, "as above", ""},
    Row{"31.09", "lombardi_silicon_holes.B", "9.93e6", "cm/s", Kind::printed, Grade::PB, "as above",
        "COMSOL's extraction gave 9.93e7 (Unit 46)"},
    Row{"31.09", "lombardi_silicon_holes.C", "8.84e5", "", Kind::printed, Grade::PB, "as above", ""},
    Row{"31.09", "lombardi_silicon_holes.lambda", "0.0317", "", Kind::printed, Grade::PB, "as above", ""},
    Row{"31.09", "lombardi_silicon_holes.delta", "2.05e14", "V/s", Kind::printed, Grade::PB, "as above", ""},
    // 31.10 energy relaxation times.
    Row{"31.10", "energy_relaxation.electron", "0.35e-12", "s", Kind::printed, Grade::PA,
        "A. Triebl, PhD thesis, TU Wien 2012, sec. 4.3.1 (\"typical values\", citing Vasicek 2009)", ""},
    Row{"31.10", "energy_relaxation.hole", "0.4e-12", "s", Kind::printed, Grade::PA,
        "Triebl 2012, sec. 4.3.1 (citing Palankovski 2000)", ""},
    Row{"31.10", "energy_balance_regularization", "1e-8", "", Kind::choice, Grade::C,
        "owner ruling 2026-10-29 (Unit 52, \"beta = 1\": regularize below a0 = 1e-8)",
        "rebound in step 4 (15.14): the row named energy_balance.canali_beta, which binds no field; the ruling fixes a0"},
    // 31.11 Scharfetter lifetimes.
    Row{"31.11", "legacy_silicon.lifetime.tau_n0", "", "s", Kind::printed, Grade::L, "not found",
        "Genius ships 1e-7; Xyce 3.95e-4 (Weaver-Nasby)"},
    Row{"31.11", "legacy_silicon.lifetime.tau_p0", "", "s", Kind::printed, Grade::L, "not found",
        "Genius 1e-7; Xyce 3.52e-5"},
    Row{"31.11", "legacy_silicon.lifetime.N_ref", "5e16", "cm^-3", Kind::printed, Grade::PC,
        "Genius Si_band_default.cc NSRHN, NSRHP (uncited)", "form cited to Fossum et al. (Gritsch sec. 4.7)"},
    // 31.13 Auger.
    Row{"31.13", "legacy_silicon.auger.Cn", "2.8e-31", "cm^6/s", Kind::printed, Grade::PB,
        "Palankovski Table 3.36; Xyce MaterialLayer.C; Genius AUGN (all uncited)",
        "the Dziewior-Schmid attribution is unconfirmed (L)"},
    Row{"31.13", "legacy_silicon.auger.Cp", "9.9e-32", "cm^6/s", Kind::printed, Grade::PB, "as above", ""},
    // 31.17 dopant levels.
    Row{"31.17", "legacy_silicon.ionization.donor_eV", "0.045", "eV", Kind::printed, Grade::PC,
        "Genius IncompleteIonization_Init (phosphorus E0)", ""},
    Row{"31.17", "legacy_silicon.ionization.acceptor_eV", "0.045", "eV", Kind::printed, Grade::PC,
        "Genius (boron E0)", ""},
    Row{"31.17", "legacy_silicon.ionization.donor_degeneracy", "2.0", "", Kind::printed, Grade::PC,
        "Genius (phosphorus GB)", ""},
    Row{"31.17", "legacy_silicon.ionization.acceptor_degeneracy", "4.0", "", Kind::printed, Grade::PC,
        "Genius (boron GB)", ""},
    // 31.18 Altermatt-Schenk, Table I (Unit 37 read; the paper's g is the reciprocal degeneracy).
    Row{"31.18", "altermatt_phosphorus.E_dop0", "45.5e-3", "eV", Kind::printed, Grade::PA,
        "Schenk, Altermatt, Schmithuesen, SISPAD 2006, Table I", ""},
    Row{"31.18", "altermatt_phosphorus.N_ref", "2.2e18", "cm^-3", Kind::printed, Grade::PA, "as above", ""},
    Row{"31.18", "altermatt_phosphorus.c", "2", "", Kind::printed, Grade::PA, "as above", ""},
    Row{"31.18", "altermatt_phosphorus.N_b", "6e18", "cm^-3", Kind::printed, Grade::PA, "as above", ""},
    Row{"31.18", "altermatt_phosphorus.d", "2.3", "", Kind::printed, Grade::PA, "as above", ""},
    Row{"31.18", "altermatt_boron.E_dop0", "44.39e-3", "eV", Kind::printed, Grade::PA, "as above", ""},
    Row{"31.18", "altermatt_boron.N_ref", "1.3e18", "cm^-3", Kind::printed, Grade::PA, "as above", ""},
    Row{"31.18", "altermatt_boron.c", "1.4", "", Kind::printed, Grade::PA, "as above", ""},
    Row{"31.18", "altermatt_boron.N_b", "4.5e18", "cm^-3", Kind::printed, Grade::PA, "as above", ""},
    Row{"31.18", "altermatt_boron.d", "2.4", "", Kind::printed, Grade::PA, "as above", ""},
    // 31.20 van Overstraeten-de Man (F4).
    Row{"31.20", "legacy_silicon.impact_ionization.electron.A", "7.03e5", "1/cm", Kind::printed, Grade::PC,
        "van Overstraeten, de Man, Solid-State Electron. 13, 583 (1970), abstract, as Allpix", ""},
    Row{"31.20", "legacy_silicon.impact_ionization.electron.B", "1.231e6", "V/cm", Kind::printed, Grade::PC,
        "as above", ""},
    Row{"31.20", "legacy_silicon.impact_ionization.hole.A_low", "1.582e6", "1/cm", Kind::printed, Grade::PC,
        "as above", ""},
    Row{"31.20", "legacy_silicon.impact_ionization.hole.B_low", "2.036e6", "V/cm", Kind::printed, Grade::PC,
        "as above", ""},
    Row{"31.20", "legacy_silicon.impact_ionization.hole.A_high", "6.71e5", "1/cm", Kind::printed, Grade::PC,
        "as above", ""},
    Row{"31.20", "legacy_silicon.impact_ionization.hole.B_high", "1.693e6", "V/cm", Kind::printed, Grade::PC,
        "as above", ""},
    Row{"31.20", "legacy_silicon.impact_ionization.hole.switch", "4.0e5", "V/cm", Kind::printed, Grade::PC,
        "as above (E0)", "the electron switch, 5e5, is immaterial: both electron branches are equal"},
    Row{"31.20", "legacy_silicon.impact_ionization.phonon_energy_eV", "0.063", "eV", Kind::printed, Grade::PB,
        "Allpix; Genius OP_PH_EN", ""},
    Row{"31.20", "impact_ionization.gamma_at_300K", "1", "", Kind::choice, Grade::C,
        "Sentaurus manual only (F4); T0 = 300 K assumed (Allpix)",
        "pins T0 = 300 K through its consequence gamma(300 K) = 1 (bit-exact: tanh(x) / tanh(x)); the form is F4's"},
    // 31.21 Kane pair (F5).
    Row{"31.21", "legacy_silicon.band_to_band.A_per_cm3_s", "", "cm^-3 s^-1", Kind::printed, Grade::L,
        "not found: none of Jiao et al. 2013 Table 1's silicon sets or Sentaurus' Hurkx defaults (4e14)", "F5"},
    Row{"31.21", "legacy_silicon.band_to_band.B_V_per_cm", "", "V/cm", Kind::printed, Grade::L,
        "not found (Sentaurus Hurkx 1.9e7; Jiao direct 1.01e8)", "F5"},
    // 31.28 thermal.
    Row{"31.28", "legacy_silicon.thermal.conductivity_W_cmK", "1.48", "W/(cm K)", Kind::printed, Grade::PC,
        "Palankovski Table 3.5 (148 W/(K m); refs Ng, Sze)", ""},
    Row{"31.28", "legacy_silicon.thermal.conductivity_exponent", "1.3", "", Kind::printed, Grade::PC,
        "Palankovski Table 3.5 (alpha = -1.3)", "shipped 1.33: a mismatch (15.10)"},
    Row{"31.28", "legacy_silicon.thermal.heat_capacity_J_cm3K", "", "J/(cm^3 K)", Kind::printed, Grade::L,
        "Ioffe NSM (cited in the header, not read)", ""},
    // 31.30 SiO2 (15.19, MI-Q26: the offset row rebound as chi_ox; the masses split; kappa and rho c added).
    Row{"31.30", "silicon_dioxide.eps_r", "", "", Kind::printed, Grade::L, "legacy EPS_OX_R; not read", ""},
    Row{"31.30", "silicon_dioxide_tunnelling.band_gap_eV", "9.0", "eV", Kind::printed, Grade::PA,
        "Gehring, PhD thesis, TU Wien 2003, Table 5.1 (from Leroy 2003)", "Unit 51"},
    Row{"31.30", "silicon_dioxide_tunnelling.electron_affinity_eV", "", "eV", Kind::printed, Grade::L,
        "chi_Si (L) less the conduction-band offset 3.0 eV of Gehring 2003, Table 5.1 (PA)",
        "1.05 = 4.05 - 3.0; graded by its L input, silicon's chi (MI-Q26)"},
    Row{"31.30", "silicon_dioxide_tunnelling.tunnel_mass_n", "0.50", "m0", Kind::choice, Grade::C,
        "owner ruling 2026-10-28 (Unit 51)", ""},
    Row{"31.30", "silicon_dioxide_tunnelling.tunnel_mass_p", "0.50", "m0", Kind::choice, Grade::C,
        "owner ruling 2026-10-28 (Unit 51)", ""},
    Row{"31.30", "silicon_dioxide.thermal.conductivity_W_cmK", "", "W/(cm K)", Kind::printed, Grade::L,
        "\"about 1.4 W/(m K) at 300 K\" (code comment; no source read)", ""},
    Row{"31.30", "silicon_dioxide.thermal.heat_capacity_J_cm3K", "", "J/(cm^3 K)", Kind::printed, Grade::L,
        "2.20 g/cm^3 x 0.74 J/(g K) (NIST-JANAF cited in the code, not read)", ""},
    // 31.27 MLDA quantization masses.
    Row{"31.27", "mlda_silicon_electron_mass", "0.98", "m0", Kind::printed, Grade::PC, "Ioffe (Unit 47)", ""},
    Row{"31.27", "mlda_silicon_hole_mass", "0.49", "m0", Kind::printed, Grade::PC, "Ioffe (Unit 47)", ""},
    // 31.31 the non-silicon sets (15.19). Ioffe = Ioffe NSM archive, Semicond/<material>/<page>.html, read
    // 2026-11-06; PA where the code cites Ioffe for the value, otherwise PC (MI-Q24); "~", "<=" and "<" qualifiers
    // in the note.
    Row{"31.31", "germanium.eps_r", "16.2", "", Kind::printed, Grade::PC, "Ioffe Ge/basic", ""},
    Row{"31.31", "germanium.Eg0_eV", "0.742", "eV", Kind::printed, Grade::PC, "Ioffe Ge/bandstr (E_g(T))", ""},
    Row{"31.31", "germanium.varshni_alpha_eV_per_K", "4.8e-4", "eV/K", Kind::printed, Grade::PC,
        "Ioffe Ge/bandstr: E_g = 0.742 - 4.8e-4 T^2/(T + 235)", ""},
    Row{"31.31", "germanium.varshni_beta_K", "235", "K", Kind::printed, Grade::PC, "Ioffe Ge/bandstr", ""},
    Row{"31.31", "germanium.Nc300", "1.0e19", "cm^-3", Kind::printed, Grade::PC, "Ioffe Ge/bandstr", ""},
    Row{"31.31", "germanium.Nv300", "5.0e18", "cm^-3", Kind::printed, Grade::PC, "Ioffe Ge/bandstr", ""},
    Row{"31.31", "germanium.electron_affinity_eV", "4.0", "eV", Kind::printed, Grade::PC, "Ioffe Ge/basic", ""},
    Row{"31.31", "germanium.electron_mobility.mu_max", "3900", "cm^2/(V s)", Kind::printed, Grade::PC,
        "Ioffe Ge/electric", "printed <= 3900"},
    Row{"31.31", "germanium.hole_mobility.mu_max", "1900", "cm^2/(V s)", Kind::printed, Grade::PC,
        "Ioffe Ge/electric", "printed <= 1900"},
    Row{"31.31", "germanium.electron_mobility.T_exponent", "-1.66", "", Kind::printed, Grade::PC,
        "Ioffe Ge/electric: mu_n ~ 4.9e7 T^-1.66, weakly doped, 77-300 K", ""},
    Row{"31.31", "germanium.hole_mobility.T_exponent", "-2.33", "", Kind::printed, Grade::PC,
        "Ioffe Ge/electric: mu_p ~ 1.05e9 T^-2.33, weakly doped, 100-300 K", ""},
    Row{"31.31", "germanium.auger.Cn", "1e-30", "cm^6/s", Kind::printed, Grade::PC, "Ioffe Ge/electric",
        "printed ~10^-30, one coefficient"},
    Row{"31.31", "germanium.auger.Cp", "1e-30", "cm^6/s", Kind::printed, Grade::PC, "Ioffe Ge/electric",
        "printed ~10^-30, one coefficient"},
    Row{"31.31", "germanium.electron_saturation.v_sat_cm_s", "", "cm/s", Kind::printed, Grade::L,
        "not printed (Ioffe Ge/electric: figures only)", ""},
    Row{"31.31", "germanium.hole_saturation.v_sat_cm_s", "", "cm/s", Kind::printed, Grade::L, "not printed", ""},
    Row{"31.31", "germanium.electron_saturation.beta", "", "", Kind::printed, Grade::L, "not printed", ""},
    Row{"31.31", "germanium.hole_saturation.beta", "", "", Kind::printed, Grade::L, "not printed", ""},
    Row{"31.31", "germanium.radiative_cm3_s", "6.41e-14", "cm^3/s", Kind::printed, Grade::PA,
        "Ioffe Ge/electric (radiative recombination coefficient at 300 K)", "the code cites Ioffe"},
    Row{"31.31", "germanium.ionization.donor_eV", "0.013", "eV", Kind::printed, Grade::PC, "Ioffe Ge/bandstr (P)",
        ""},
    Row{"31.31", "germanium.ionization.acceptor_eV", "0.011", "eV", Kind::printed, Grade::PC,
        "Ioffe Ge/bandstr (B)", ""},
    Row{"31.31", "gallium_arsenide.eps_r", "12.9", "", Kind::printed, Grade::PC, "Ioffe GaAs/basic", ""},
    Row{"31.31", "gallium_arsenide.Eg0_eV", "1.519", "eV", Kind::printed, Grade::PC, "Ioffe GaAs/bandstr", ""},
    Row{"31.31", "gallium_arsenide.varshni_alpha_eV_per_K", "5.405e-4", "eV/K", Kind::printed, Grade::PC,
        "Ioffe GaAs/bandstr: E_g = 1.519 - 5.405e-4 T^2/(T + 204)", ""},
    Row{"31.31", "gallium_arsenide.varshni_beta_K", "204", "K", Kind::printed, Grade::PC, "Ioffe GaAs/bandstr",
        ""},
    Row{"31.31", "gallium_arsenide.Nc300", "4.7e17", "cm^-3", Kind::printed, Grade::PC, "Ioffe GaAs/bandstr", ""},
    Row{"31.31", "gallium_arsenide.Nv300", "9.0e18", "cm^-3", Kind::printed, Grade::PC, "Ioffe GaAs/bandstr", ""},
    Row{"31.31", "gallium_arsenide.electron_affinity_eV", "4.07", "eV", Kind::printed, Grade::PC,
        "Ioffe GaAs/basic", ""},
    Row{"31.31", "gallium_arsenide.electron_mobility.mu_max", "8500", "cm^2/(V s)", Kind::printed, Grade::PC,
        "Ioffe GaAs/electric", "printed <= 8500"},
    Row{"31.31", "gallium_arsenide.hole_mobility.mu_max", "400", "cm^2/(V s)", Kind::printed, Grade::PC,
        "Ioffe GaAs/electric", "printed <= 400"},
    Row{"31.31", "gallium_arsenide.auger.Cn", "1e-30", "cm^6/s", Kind::printed, Grade::PC, "Ioffe GaAs/electric",
        "printed ~10^-30 at 300 K, one coefficient"},
    Row{"31.31", "gallium_arsenide.auger.Cp", "1e-30", "cm^6/s", Kind::printed, Grade::PC, "Ioffe GaAs/electric",
        "printed ~10^-30 at 300 K, one coefficient"},
    Row{"31.31", "gallium_arsenide.electron_saturation.v_sat_cm_s", "", "cm/s", Kind::printed, Grade::L,
        "not printed (Ioffe GaAs/electric: figures only)", ""},
    Row{"31.31", "gallium_arsenide.hole_saturation.v_sat_cm_s", "", "cm/s", Kind::printed, Grade::L,
        "not printed", ""},
    Row{"31.31", "gallium_arsenide.electron_saturation.beta", "", "", Kind::printed, Grade::L, "not printed", ""},
    Row{"31.31", "gallium_arsenide.hole_saturation.beta", "", "", Kind::printed, Grade::L, "not printed", ""},
    Row{"31.31", "gallium_arsenide.radiative_cm3_s", "7.2e-10", "cm^3/s", Kind::printed, Grade::PC,
        "Ioffe GaAs/electric (Varshni 1967, 300 K)", "the code cites nextnano (not read)"},
    Row{"31.31", "gallium_arsenide.ionization.donor_eV", "0.006", "eV", Kind::printed, Grade::PC,
        "Ioffe GaAs/bandstr (Si; Milnes 1973)", "printed ~0.006; the code cites Phys. Rev. B 43, 14734 (not read)"},
    Row{"31.31", "gallium_arsenide.ionization.acceptor_eV", "0.02", "eV", Kind::printed, Grade::PC,
        "Ioffe GaAs/bandstr (C; Milnes 1973)", "printed ~0.02; the code cites Phys. Rev. B 43, 14734 (not read)"},
    Row{"31.31", "indium_gallium_arsenide.eps_r", "13.9", "", Kind::printed, Grade::PC,
        "Ioffe GaInAs/basic (Ga0.47In0.53As)", ""},
    Row{"31.31", "indium_gallium_arsenide.Eg0_eV", "", "eV", Kind::printed, Grade::L,
        "Ioffe GaInAs/basic prints 0.74 eV at 300 K and Goetz's 2 K law 0.4105 + 0.6337x + 0.475x^2", ""},
    Row{"31.31", "indium_gallium_arsenide.varshni_alpha_eV_per_K", "", "eV/K", Kind::printed, Grade::L,
        "Ioffe GaInAs/bandstr prints Paul et al.'s E_g(x, T), not a Varshni pair", ""},
    Row{"31.31", "indium_gallium_arsenide.varshni_beta_K", "", "K", Kind::printed, Grade::L, "as above", ""},
    Row{"31.31", "indium_gallium_arsenide.Nc300", "2.1e17", "cm^-3", Kind::printed, Grade::PC,
        "Ioffe GaInAs/bandstr", ""},
    Row{"31.31", "indium_gallium_arsenide.Nv300", "7.7e18", "cm^-3", Kind::printed, Grade::PC,
        "Ioffe GaInAs/bandstr", ""},
    Row{"31.31", "indium_gallium_arsenide.electron_affinity_eV", "4.5", "eV", Kind::printed, Grade::PC,
        "Ioffe GaInAs/basic", ""},
    Row{"31.31", "indium_gallium_arsenide.electron_mobility.mu_max", "12e3", "cm^2/(V s)", Kind::printed,
        Grade::PC, "Ioffe GaInAs/electric", "printed < 12 10^3"},
    Row{"31.31", "indium_gallium_arsenide.hole_mobility.mu_max", "300", "cm^2/(V s)", Kind::printed, Grade::PC,
        "Ioffe GaInAs/electric", "printed < 300"},
    Row{"31.31", "indium_gallium_arsenide.auger.Cn", "7e-29", "cm^6/s", Kind::printed, Grade::PC,
        "Ioffe GaInAs/electric", "one coefficient"},
    Row{"31.31", "indium_gallium_arsenide.auger.Cp", "7e-29", "cm^6/s", Kind::printed, Grade::PC,
        "Ioffe GaInAs/electric", "one coefficient"},
    Row{"31.31", "indium_gallium_arsenide.electron_saturation.v_sat_cm_s", "", "cm/s", Kind::printed, Grade::L,
        "not printed", ""},
    Row{"31.31", "indium_gallium_arsenide.hole_saturation.v_sat_cm_s", "", "cm/s", Kind::printed, Grade::L,
        "not printed", ""},
    Row{"31.31", "indium_gallium_arsenide.electron_saturation.beta", "", "", Kind::printed, Grade::L,
        "not printed", ""},
    Row{"31.31", "indium_gallium_arsenide.hole_saturation.beta", "", "", Kind::printed, Grade::L, "not printed",
        ""},
    Row{"31.31", "indium_gallium_arsenide.radiative_cm3_s", "0.96e-10", "cm^3/s", Kind::printed, Grade::PA,
        "Ioffe GaInAs/electric", "the code cites Ioffe"},
    Row{"31.31", "silicon_carbide_4h.eps_r", "", "", Kind::printed, Grade::L,
        "Ioffe SiC/basic: for 4H the value of 6H-SiC is usually used (9.66 ordinary, 10.03 extraordinary)", ""},
    Row{"31.31", "silicon_carbide_4h.Eg0_eV", "", "eV", Kind::printed, Grade::L,
        "Ioffe SiC/bandstr prints E_g 3.23 eV at 300 K and E_g(0) - 6.5e-4 T^2/(T + 1300), E_g(0) not printed", ""},
    Row{"31.31", "silicon_carbide_4h.varshni_alpha_eV_per_K", "6.5e-4", "eV/K", Kind::printed, Grade::PC,
        "Ioffe SiC/bandstr (Goldberg et al. 2001)", ""},
    Row{"31.31", "silicon_carbide_4h.varshni_beta_K", "1300", "K", Kind::printed, Grade::PC,
        "Ioffe SiC/bandstr (Goldberg et al. 2001)", ""},
    Row{"31.31", "silicon_carbide_4h.Nc300", "", "cm^-3", Kind::printed, Grade::L,
        "Ioffe SiC/bandstr prints N_c ~ 3.25e15 T^(3/2) (1.69e19 at 300 K)", ""},
    Row{"31.31", "silicon_carbide_4h.Nv300", "", "cm^-3", Kind::printed, Grade::L,
        "Ioffe SiC/bandstr prints N_v ~ 4.8e15 T^(3/2) (2.49e19 at 300 K)", ""},
    Row{"31.31", "silicon_carbide_4h.electron_affinity_eV", "", "eV", Kind::printed, Grade::L,
        "not printed (the legacy cites Baliga)", ""},
    Row{"31.31", "silicon_carbide_4h.electron_mobility.mu_max", "900", "cm^2/(V s)", Kind::printed, Grade::PC,
        "Ioffe SiC/ebasic", "printed <= 900"},
    Row{"31.31", "silicon_carbide_4h.hole_mobility.mu_max", "120", "cm^2/(V s)", Kind::printed, Grade::PC,
        "Ioffe SiC/ebasic", "printed <= 120"},
    Row{"31.31", "silicon_carbide_4h.electron_mobility.mu_min", "", "cm^2/(V s)", Kind::printed, Grade::L,
        "not printed (Caughey-Thomas fit carried by the legacy)", ""},
    Row{"31.31", "silicon_carbide_4h.hole_mobility.mu_min", "", "cm^2/(V s)", Kind::printed, Grade::L,
        "not printed", ""},
    Row{"31.31", "silicon_carbide_4h.electron_mobility.N_ref", "", "cm^-3", Kind::printed, Grade::L,
        "not printed", ""},
    Row{"31.31", "silicon_carbide_4h.hole_mobility.N_ref", "", "cm^-3", Kind::printed, Grade::L, "not printed", ""},
    Row{"31.31", "silicon_carbide_4h.electron_mobility.alpha", "", "", Kind::printed, Grade::L, "not printed", ""},
    Row{"31.31", "silicon_carbide_4h.hole_mobility.alpha", "", "", Kind::printed, Grade::L, "not printed", ""},
    Row{"31.31", "silicon_carbide_4h.electron_mobility.T_exponent", "", "", Kind::printed, Grade::L,
        "not printed", ""},
    Row{"31.31", "silicon_carbide_4h.hole_mobility.T_exponent", "", "", Kind::printed, Grade::L, "not printed",
        ""},
    Row{"31.31", "silicon_carbide_4h.lifetime.tau_n0", "1e-9", "s", Kind::printed, Grade::PC,
        "Ioffe SiC/recombination: longest electron lifetime, pure p-type, 300 K", "printed ~10^-9"},
    Row{"31.31", "silicon_carbide_4h.lifetime.tau_p0", "6.0e-7", "s", Kind::printed, Grade::PC,
        "Ioffe SiC/recombination: longest hole lifetime, pure n-type, 300 K", "printed ~6.0e-7"},
    Row{"31.31", "silicon_carbide_4h.lifetime.N_ref", "", "cm^-3", Kind::printed, Grade::L, "not printed", ""},
    Row{"31.31", "silicon_carbide_4h.auger.Cn", "5.0e-31", "cm^6/s", Kind::printed, Grade::PC,
        "Ioffe SiC/recombination (Galeskas et al. 1997)", ""},
    Row{"31.31", "silicon_carbide_4h.auger.Cp", "2.0e-31", "cm^6/s", Kind::printed, Grade::PC,
        "Ioffe SiC/recombination (Galeskas et al. 1997)", ""},
    Row{"31.31", "silicon_carbide_4h.electron_saturation.v_sat_cm_s", "", "cm/s", Kind::printed, Grade::L,
        "not printed", ""},
    Row{"31.31", "silicon_carbide_4h.hole_saturation.v_sat_cm_s", "", "cm/s", Kind::printed, Grade::L,
        "not printed", ""},
    Row{"31.31", "silicon_carbide_4h.electron_saturation.beta", "", "", Kind::printed, Grade::L, "not printed", ""},
    Row{"31.31", "silicon_carbide_4h.hole_saturation.beta", "", "", Kind::printed, Grade::L, "not printed", ""},
    Row{"31.31", "silicon_carbide_4h.radiative_cm3_s", "1.5e-12", "cm^3/s", Kind::printed, Grade::PA,
        "Ioffe SiC/basic and SiC/recombination (Galeskas et al. 1997, estimation)", "the code cites Ioffe"},
    Row{"31.31", "silicon_carbide_4h.ionization.donor_eV", "", "eV", Kind::printed, Grade::L,
        "Ayalew, TU Wien thesis (cited in the code, not read)", ""},
    Row{"31.31", "silicon_carbide_4h.ionization.acceptor_eV", "", "eV", Kind::printed, Grade::L,
        "Ayalew, TU Wien thesis (cited in the code, not read)", ""},
    Row{"31.31", "silicon_carbide_4h.thermal.conductivity_W_cmK", "3.7", "W/(cm K)", Kind::printed, Grade::PC,
        "Ioffe SiC/thermal (4H, Goldberg et al. 2001)", "the code cites the legacy materials.py"},
    Row{"31.31", "silicon_carbide_4h.thermal.heat_capacity_J_cm3K", "", "J/(cm^3 K)", Kind::printed, Grade::L,
        "Ioffe density 3.211 g/cm^3 (4H, SiC/basic) and specific heat 0.69 J/(g K) (6H, usually used; SiC/thermal)",
        "their product is 2.21559 (not printed); shipped 2.2156"},
    Row{"31.31", "algaas.radiative_cm3_s", "1.8e-10", "cm^3/s", Kind::printed, Grade::PA, "Ioffe AlGaAs/recomb",
        "printed ~1.8e-10 at 300 K; the code cites Ioffe"},
};

// The printed temperature laws behind the law_300K rows: value(T) = coefficient T^exponent, typed from the source
// (Allpix Squared 3.3 manual's Canali section, Jacoboni et al. 1977 Table 5). The table check judges the shipped
// value against the law at 300 K to the law's printed digits: the sum of the value's changes for half a unit in
// the last printed digit of each (15.9 "Parameter register" with the law_300K kind above).
struct Law {
    std::string_view binding;
    std::string_view coefficient;
    std::string_view exponent;
};
inline constexpr std::array laws{
    Law{"legacy_silicon.electron_saturation.v_sat_cm_s", "1.53e9", "-0.87"},
    Law{"legacy_silicon.hole_saturation.v_sat_cm_s", "1.62e8", "-0.52"},
    Law{"legacy_silicon.electron_saturation.beta", "2.57e-2", "0.66"},
    Law{"legacy_silicon.hole_saturation.beta", "0.46", "0.17"},
};

// The Al_x Ga_1-x As composition laws of 31.31 (15.19, MI-Q23), typed from Ioffe AlGaAs/basic, bandstr and ebasic
// (300 K): value(x) = c0 + c1 x + c2 x^2 (polynomial), or scale (c0 + c1 x)^(3/2) (dos_power). The check judges
// algaas_parameters(x) to the change of the law for half a unit in the last printed digit of each coefficient
// (summed), on 0 <= x <= x_max. Grade PC.
enum class LawForm { polynomial, dos_power };
struct CompositionLaw {
    std::string_view binding;  // "algaas.<field>", or "algaas.band_gap_300K" for E_g(300 K) of the set
    LawForm form;
    std::string_view scale;    // dos_power only
    std::array<std::string_view, 3> c;  // empty: absent
    std::string_view x_max;
    std::string_view unit;
    std::string_view source;
};
inline constexpr std::array composition_laws{
    CompositionLaw{"algaas.eps_r", LawForm::polynomial, "", {"12.90", "-2.84", ""}, "0.45", "",
                   "Ioffe AlGaAs/basic: 12.90-2.84x"},
    CompositionLaw{"algaas.band_gap_300K", LawForm::polynomial, "", {"1.424", "1.247", ""}, "0.45", "eV",
                   "Ioffe AlGaAs/bandstr: x<0.45 1.424+1.247x eV"},
    CompositionLaw{"algaas.electron_affinity_eV", LawForm::polynomial, "", {"4.07", "-1.1", ""}, "0.45", "eV",
                   "Ioffe AlGaAs/basic: 4.07-1.1x eV (x<0.45)"},
    CompositionLaw{"algaas.Nc300", LawForm::dos_power, "2.5e19", {"0.063", "0.083", ""}, "0.41", "cm^-3",
                   "Ioffe AlGaAs/bandstr: x<0.41 2.5e19 (0.063+0.083x)^(3/2)"},
    CompositionLaw{"algaas.Nv300", LawForm::dos_power, "2.5e19", {"0.51", "0.25", ""}, "1", "cm^-3",
                   "Ioffe AlGaAs/bandstr: 2.5e19 (0.51+0.25x)^(3/2)"},
    CompositionLaw{"algaas.electron_mobility.mu_max", LawForm::polynomial, "", {"8e3", "-2.2e4", "1e4"}, "0.45",
                   "cm^2/(V s)", "Ioffe AlGaAs/ebasic: 0<x<0.45 8e3-2.2e4x+1e4x^2"},
    CompositionLaw{"algaas.hole_mobility.mu_max", LawForm::polynomial, "", {"370", "-970", "740"}, "1",
                   "cm^2/(V s)", "Ioffe AlGaAs/ebasic: 370-970x+740x^2"},
};

}  // namespace NiTCAD::verify::reference::model_parameters
