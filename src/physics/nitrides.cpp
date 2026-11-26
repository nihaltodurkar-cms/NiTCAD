#include "NiTCAD/physics/nitrides.hpp"

#include <cmath>
#include <format>
#include <limits>

#include "temperature.hpp"

namespace NiTCAD::physics {
namespace {

// Nc or Nv at 300 K from a density-of-states mass [m0], the engine's own constant.
double dos_300(double mass) { return detail::dos_constant * std::pow(mass, 1.5) * std::pow(300.0, 1.5); }

double varshni_gap_300(double Eg0, double alpha, double beta) { return Eg0 - alpha * 300.0 * 300.0 / (300.0 + beta); }

// A lifetime with no doping dependence: tau0 / (1 + N / N_ref) is tau0 to the last bit.
constexpr double constant_lifetime_N_ref = std::numeric_limits<double>::max();

// Constant hole mobility (amendment 6): mu_min = mu_max.
constexpr CaugheyThomasParameters constant_mobility(double mu) {
    return {.mu_min = mu, .mu_max = mu, .N_ref = 1.0e17, .alpha = 1.0, .T_exponent = 0.0};
}

SemiconductorParameters nitride(double eps, double Eg0, double alpha, double beta, double me, double mh,
                                CaugheyThomasParameters electron, double mu_h, double vsat_n, double beta_n,
                                double tau_n, double tau_p, double Cn, double Cp, double B, double donor,
                                double acceptor, double k300, double delta, double heat_J_cm3K) {
    SemiconductorParameters p{
        .eps_r = eps,
        .Eg0_eV = Eg0,
        .varshni_alpha_eV_per_K = alpha,
        .varshni_beta_K = beta,
        .Nc300 = dos_300(me),
        .Nv300 = dos_300(mh),
        .electron_affinity_eV = gan_affinity_eV,
        .electron_mobility = electron,
        .hole_mobility = constant_mobility(mu_h),
        .lifetime = {.tau_n0 = tau_n, .tau_p0 = tau_p, .N_ref = constant_lifetime_N_ref},
        .auger = {.Cn = Cn, .Cp = Cp},
        .bandgap_narrowing = no_bandgap_narrowing,  // off (amendment 6)
        .electron_saturation = {.v_sat_cm_s = vsat_n, .beta = beta_n},
        .hole_saturation = {.v_sat_cm_s = 7.0e6, .beta = 1.0},
        .radiative_cm3_s = B,
        .ionization = {.donor_eV = donor, .acceptor_eV = acceptor, .donor_degeneracy = 2.0,
                       .acceptor_degeneracy = 4.0},
        .richardson = {.electron = 0.0, .hole = 0.0},
    };
    // NiTCAD's kappa(T) = kappa300 (T / 300)^(-conductivity_exponent); the owner's delta is the
    // exponent itself.
    p.thermal = {.conductivity_W_cmK = k300, .conductivity_exponent = -delta, .heat_capacity_J_cm3K = heat_J_cm3K};
    return p;
}

// The 300 K gap and the Q_c alignment against GaN.
double aligned_affinity(const SemiconductorParameters& p, double share) {
    const SemiconductorParameters g = gan_parameters();
    const double gap_gan = varshni_gap_300(g.Eg0_eV, g.varshni_alpha_eV_per_K, g.varshni_beta_K);
    const double gap = varshni_gap_300(p.Eg0_eV, p.varshni_alpha_eV_per_K, p.varshni_beta_K);
    return gan_affinity_eV - share * (gap - gap_gan);
}

}  // namespace

SemiconductorParameters gan_parameters() {
    SemiconductorParameters p = nitride(
        8.9, 3.47, 7.7e-4, 600.0, 0.20, 1.5,
        {.mu_min = 80.0, .mu_max = 1405.0, .N_ref = 7.78e16, .alpha = 0.71, .T_exponent = -2.85,
         .mu_min_T_exponent = -0.20, .N_ref_T_exponent = 1.30, .alpha_T_exponent = 0.31},
        170.0, 2.5e7, 1.0, 10e-9, 20e-9, 1.25e-31, 5.0e-32, 1.1e-8, 0.017, 0.17, 1.3, -0.6, 0.49 * 6.15);
    // Measured Chynoweth coefficients (owner, amendment 6), one law over all fields (low and high
    // the same, so the switch field, which must be positive, changes nothing); no phonon
    // temperature factor given (0: none).
    p.impact_ionization = {.electron = {.A_low_per_cm = 4.48e8, .B_low_V_per_cm = 3.39e7, .A_high_per_cm = 4.48e8,
                                        .B_high_V_per_cm = 3.39e7, .switch_V_per_cm = 1.0},
                           .hole = {.A_low_per_cm = 7.13e6, .B_low_V_per_cm = 1.46e7, .A_high_per_cm = 7.13e6,
                                    .B_high_V_per_cm = 1.46e7, .switch_V_per_cm = 1.0},
                           .phonon_energy_eV = 0.0};
    return p;
}

SemiconductorParameters aln_parameters() {
    SemiconductorParameters p = nitride(
        8.5, 6.10, 2.63e-3, 2082.0, 0.40, 2.7,
        // alpha +0.18: the owner's provisional sign correction of the source's -0.18.
        {.mu_min = 29.0, .mu_max = 683.0, .N_ref = 5.0e17, .alpha = 0.18, .T_exponent = 0.80,
         .mu_min_T_exponent = -3.21, .N_ref_T_exponent = 1.21, .alpha_T_exponent = 0.31},
        // beta_n 1.0, not the source's 0.45: the engine needs beta >= 1 (owner ruling 2026-11-27).
        14.0, 3.7e7, 1.0, 5e-9, 10e-9, 5.0e-32, 5.0e-32, 0.4e-10, 0.25, 0.51, 2.85, -1.2, 0.6 * 3.23);
    p.electron_affinity_eV = aligned_affinity(p, algan_conduction_share);
    return p;
}

SemiconductorParameters inn_parameters() {
    SemiconductorParameters p = nitride(
        15.3, 0.69, 4.14e-4, 154.0, 0.11, 1.65,
        {.mu_min = 500.0, .mu_max = 10400.0, .N_ref = 3.4e17, .alpha = 0.70, .T_exponent = -3.70,
         .mu_min_T_exponent = 2.39, .N_ref_T_exponent = -0.33, .alpha_T_exponent = 0.31},
        // Donor 0: fully ionized, the owner's fallback (no InN:Si level fixed), not a measured 0 eV.
        220.0, 5.0e7, 1.0, 10e-9, 20e-9, 7.0e-30, 6.0e-30, 2.0e-10, 0.0, 0.06, 0.45, 0.0, 0.32 * 6.81);
    p.electron_affinity_eV = aligned_affinity(p, ingan_conduction_share);
    return p;
}

namespace {

double lerp(double a, double b, double x) { return (1.0 - x) * a + x * b; }

// The alloy between GaN and binary b (header). x in (0, 1).
SemiconductorParameters nitride_alloy(const SemiconductorParameters& b, double x, double gap_bowing,
                                      double mass_bowing, double share, double C_mu) {
    const SemiconductorParameters a = gan_parameters();
    const double gap_a = varshni_gap_300(a.Eg0_eV, a.varshni_alpha_eV_per_K, a.varshni_beta_K);
    const double gap_b = varshni_gap_300(b.Eg0_eV, b.varshni_alpha_eV_per_K, b.varshni_beta_K);
    const double gap = lerp(gap_a, gap_b, x) - gap_bowing * x * (1.0 - x);
    // GaN's masses back from its Nc and Nv (the engine's constant both ways).
    const auto mass = [](double N300) { return std::pow(N300 / dos_300(1.0), 2.0 / 3.0); };
    const double me = lerp(mass(a.Nc300), mass(b.Nc300), x) - mass_bowing * x * (1.0 - x);

    SemiconductorParameters p = a;
    p.eps_r = lerp(a.eps_r, b.eps_r, x);
    // The set's own Varshni is GaN's; Eg0 puts the 300 K gap at the bowed value, and alloy_gap adds
    // x (d_b - d_a), so the gap's change from 300 K is linear in x.
    p.Eg0_eV = gap + a.varshni_alpha_eV_per_K * 300.0 * 300.0 / (300.0 + a.varshni_beta_K);
    p.alloy_gap = {.fraction = x, .alpha_a_eV_per_K = a.varshni_alpha_eV_per_K, .beta_a_K = a.varshni_beta_K,
                   .alpha_b_eV_per_K = b.varshni_alpha_eV_per_K, .beta_b_K = b.varshni_beta_K};
    p.Nc300 = dos_300(me);
    p.Nv300 = lerp(a.Nv300, b.Nv300, x);  // m_h^(3/2) linear
    p.electron_affinity_eV = gan_affinity_eV - share * (gap - gap_a);
    p.alloy_mobility = {.fraction = x, .C_mu_electron = C_mu, .electron_b = b.electron_mobility,
                        .hole_b = b.hole_mobility};
    p.lifetime.tau_n0 = lerp(a.lifetime.tau_n0, b.lifetime.tau_n0, x);
    p.lifetime.tau_p0 = lerp(a.lifetime.tau_p0, b.lifetime.tau_p0, x);
    p.auger = {.Cn = lerp(a.auger.Cn, b.auger.Cn, x), .Cp = lerp(a.auger.Cp, b.auger.Cp, x)};
    p.radiative_cm3_s = lerp(a.radiative_cm3_s, b.radiative_cm3_s, x);
    p.electron_saturation = {.v_sat_cm_s = lerp(a.electron_saturation.v_sat_cm_s, b.electron_saturation.v_sat_cm_s, x),
                             .beta = lerp(a.electron_saturation.beta, b.electron_saturation.beta, x)};
    p.ionization.donor_eV = lerp(a.ionization.donor_eV, b.ionization.donor_eV, x);
    p.ionization.acceptor_eV = lerp(a.ionization.acceptor_eV, b.ionization.acceptor_eV, x);
    p.impact_ionization = b.impact_ionization;  // off: AlN and InN carry none (owner ruling 5)
    // Thermal conductivity off (owner ruling 4): alloy scattering makes a linear rule wrong.
    p.thermal.conductivity_W_cmK = 0.0;
    p.thermal.conductivity_exponent = 0.0;
    p.thermal.heat_capacity_J_cm3K = lerp(a.thermal.heat_capacity_J_cm3K, b.thermal.heat_capacity_J_cm3K, x);
    return p;
}

std::expected<SemiconductorParameters, base::Error> alloy(const char* name, double x,
                                                          SemiconductorParameters (*end_b)(), double gap_bowing,
                                                          double mass_bowing, double share, double C_mu) {
    if (!(std::isfinite(x) && x >= 0.0 && x <= 1.0)) {
        return std::unexpected(base::Error{base::ErrorCode::invalid_input,
                                           std::format("the {} fraction {} lies outside [0, 1]", name, x)});
    }
    if (x == 0.0) return gan_parameters();
    if (x == 1.0) return end_b();
    return nitride_alloy(end_b(), x, gap_bowing, mass_bowing, share, C_mu);
}

}  // namespace

std::expected<SemiconductorParameters, base::Error> algan_parameters(double x) {
    return alloy("AlGaN Al", x, aln_parameters, 0.8, 4.8e-3, algan_conduction_share, 40.0);
}

std::expected<SemiconductorParameters, base::Error> ingan_parameters(double x) {
    return alloy("InGaN In", x, inn_parameters, 1.4, 0.0, ingan_conduction_share, 1.0e6);
}

}  // namespace NiTCAD::physics
