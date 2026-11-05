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
    Row{"31.08", "legacy_silicon.electron_saturation.v_sat_cm_s", "1.0704e7", "cm/s", Kind::law_300K, Grade::PB,
        "Jacoboni et al., Solid-State Electron. 20, 77 (1977), Table 5, v_m = 1.53e9 T^-0.87, as Allpix",
        "1.53e9 * 300^-0.87 = 1.0704e7; Canali's own eq 2a, 1.43e9 T^-0.87, gives 1.0005e7 (F2)"},
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
    Row{"31.10", "energy_balance.canali_beta", "1", "", Kind::choice, Grade::C, "owner ruling 2026-10-29 (Unit 52)",
        ""},
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
    Row{"31.20", "impact_ionization.temperature_factor_form", "", "", Kind::choice, Grade::C,
        "Sentaurus manual only (F4); T0 = 300 K assumed (Allpix)", ""},
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
    // 31.30 SiO2.
    Row{"31.30", "silicon_dioxide.eps_r", "", "", Kind::printed, Grade::L, "legacy EPS_OX_R; not read", ""},
    Row{"31.30", "silicon_dioxide_tunnelling.gap_eV", "9.0", "eV", Kind::printed, Grade::PA,
        "Gehring, PhD thesis, TU Wien 2003, Table 5.1 (from Leroy 2003)", "Unit 51"},
    Row{"31.30", "silicon_dioxide_tunnelling.conduction_offset_eV", "3.0", "eV", Kind::printed, Grade::PA,
        "Gehring 2003, Table 5.1", "Unit 51"},
    Row{"31.30", "silicon_dioxide_tunnelling.masses_m0", "0.50", "m0", Kind::choice, Grade::C,
        "owner ruling 2026-10-28 (Unit 51)", ""},
    // 31.27 MLDA quantization masses.
    Row{"31.27", "mlda_silicon_electron_mass", "0.98", "m0", Kind::printed, Grade::PC, "Ioffe (Unit 47)", ""},
    Row{"31.27", "mlda_silicon_hole_mass", "0.49", "m0", Kind::printed, Grade::PC, "Ioffe (Unit 47)", ""},
};

}  // namespace NiTCAD::verify::reference::model_parameters
