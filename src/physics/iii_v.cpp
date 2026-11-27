#include "NiTCAD/physics/iii_v.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>

#include "temperature.hpp"

namespace NiTCAD::physics {
namespace {

double dos_300(double mass) { return detail::dos_constant * std::pow(mass, 1.5) * std::pow(300.0, 1.5); }
double varshni_shift_300(double alpha, double beta) { return alpha * 300.0 * 300.0 / (300.0 + beta); }
double lerp(double a, double b, double x) { return (1.0 - x) * a + x * b; }
double hole_dos_mass(double m_hh, double m_lh) { return std::pow(std::pow(m_hh, 1.5) + std::pow(m_lh, 1.5), 2.0 / 3.0); }

constexpr double constant_lifetime_N_ref = std::numeric_limits<double>::max();  // tau0 to the last bit

// Sotoodeh's form: mu_max (300 / T)^theta1, N_ref (T / 300)^theta2, mu_min and alpha constant.
constexpr CaugheyThomasParameters sotoodeh(double mu_max, double mu_min, double N_ref, double alpha, double theta1,
                                           double theta2) {
    return {.mu_min = mu_min, .mu_max = mu_max, .N_ref = N_ref, .alpha = alpha, .T_exponent = -theta1,
            .N_ref_T_exponent = theta2};
}

// One binary as the owner gave it (amendment 7).
struct Binary {
    double eps;
    double Eg0, alpha, beta;  // the fundamental valley's gap at 0 K and Varshni [eV, eV/K, K]
    double satellite_Eg0;     // the other valley (Gamma for GaP and AlAs, X for InAs) at 0 K; 0: none
    double m_e;               // fundamental valley's density-of-states mass [m0]
    double satellite_m_e;     // the satellite valley's [m0]; 0: none
    double gamma1, gamma2;    // Luttinger
    double vbo;               // valence-band offset [eV], M1's scale (InSb 0)
    CaugheyThomasParameters electron, hole;
    double vsat_n, A_n, vsat_p, A_p;
    double tau_n, tau_p, Cn, Cp, B;
    double k300_W_cmK, k_exponent, cp_J_gK, rho_g_cm3;
    bool direct;
    bool tunnelling;  // ships the direct-gap tunnelling masses (MAD10: gaas, ingaas, inp, inas)
    double m_hh() const { return 1.0 / (gamma1 - 2.0 * gamma2); }
    double m_lh() const { return 1.0 / (gamma1 + 2.0 * gamma2); }
    double gap_300() const { return Eg0 - varshni_shift_300(alpha, beta); }
};

// X ellipsoids: three valleys, (m_l m_t^2)^(1/3) 3^(2/3).
double x_valleys(double m_l, double m_t) { return std::cbrt(m_l * m_t * m_t) * std::pow(3.0, 2.0 / 3.0); }

const Binary gaas_2001_set{12.9, 1.519, 5.405e-4, 204.0, 0.0, 0.067, 0.0, 6.98, 2.06, -0.80,
                       sotoodeh(9400.0, 500.0, 6.0e16, 0.394, 2.10, 3.0), sotoodeh(491.5, 20.0, 1.48e17, 0.38, 2.20, 3.0),
                       7.2e6, 0.44, 9.0e6, 0.59, 1.8e-7, 3.0e-8, 1.0e-30, 1.0e-30, 7.2e-10,
                       0.55, 1.25, 0.327, 5.3174, true, true};
const Binary inp_set{12.50, 1.4236, 3.63e-4, 162.0, 0.0, 0.0795, 0.0, 5.08, 1.60, -0.94,
                 sotoodeh(5200.0, 400.0, 3.0e17, 0.47, 2.00, 3.25), sotoodeh(170.0, 10.0, 4.87e17, 0.62, 2.00, 3.0),
                 6.7e6, 0.68, 5.0e6, 0.30, 2.0e-7, 2.0e-7, 9.0e-31, 9.0e-31, 6.25e-10,
                 0.68, 1.2, 0.31, 4.79, true, true};
// InAs's X valley (for InAlAs): 1.433 eV at 0 K, 0.64 m0 for the three valleys together (owner H).
const Binary inas_set{15.15, 0.417, 2.76e-4, 93.0, 1.433, 0.026, 0.64, 20.0, 8.5, -0.59,
                  sotoodeh(34000.0, 1000.0, 1.1e18, 0.32, 1.57, 3.0), sotoodeh(530.0, 20.0, 1.1e17, 0.46, 2.30, 3.0),
                  9.0e6, 0.57, 5.0e6, 0.30, 3.25e-7, 3.25e-7, 1.2e-26, 1.2e-26, 1.1e-10,
                  0.27, 1.2, 0.25, 5.67, true, true};
const Binary gap_set{11.10, 2.350, 5.771e-4, 372.0, 2.886, x_valleys(2.00, 0.253), 0.13, 4.05, 0.49, -1.27,
                 sotoodeh(152.0, 10.0, 4.4e18, 0.80, 1.60, 0.71), sotoodeh(147.0, 10.0, 1.0e18, 0.85, 1.98, 0.0),
                 8.8e6, 0.30, 5.0e6, 0.30, 1.0e-7, 1.0e-7, 0.0, 0.0, 0.0,
                 0.55, 1.2, 0.33, 4.13, false, false};
const Binary alas_set{10.06, 2.240, 7.0e-4, 530.0, 3.099, x_valleys(0.97, 0.22), 0.15, 3.76, 0.82, -1.33,
                  sotoodeh(400.0, 10.0, 5.46e17, 1.00, 2.10, 3.0), sotoodeh(200.0, 10.0, 3.84e17, 0.488, 2.24, 3.0),
                  8.5e6, 0.55, 8.0e6, 0.30, 1.0e-7, 1.0e-7, 0.0, 0.0, 0.0,
                  0.91, 1.2, 0.50, 3.73, false, false};
const Binary gasb_set{15.70, 0.812, 4.17e-4, 140.0, 0.0, 0.039, 0.0, 13.4, 4.7, -0.03,
                  sotoodeh(5650.0, 1050.0, 2.8e17, 1.05, 2.00, 2.8), sotoodeh(875.0, 190.0, 9.0e17, 0.65, 1.70, 2.7),
                  6.0e6, 0.0, 3.0e6, 0.0, 1.0e-7, 1.0e-7, 5.0e-30, 5.0e-30, 1.0e-10,
                  0.32, 1.6, 0.25, 5.61, true, false};
// InSb's mobility exponents provisional (owner); SRH off (its 300 K lifetime is Auger-limited).
const Binary insb_set{17.70, 0.235, 3.2e-4, 170.0, 0.0, 0.0135, 0.0, 34.8, 15.5, 0.0,
                  sotoodeh(78000.0, 5000.0, 7.0e16, 0.70, 1.66, 3.0), sotoodeh(750.0, 100.0, 6.0e17, 0.60, 1.80, 3.0),
                  5.5e7, 0.0, 7.0e6, 0.0, 0.0, 0.0, 5.0e-26, 5.0e-26, 5.0e-11,
                  0.18, 1.2, 0.20, 5.78, true, false};

double affinity(double vbo, double gap300) {
    return iii_v_affinity_anchor_eV + gaas_2001_set.gap_300() + gaas_2001_set.vbo - vbo - gap300;
}

SemiconductorParameters make(const Binary& b) {
    SemiconductorParameters p{
        .eps_r = b.eps,
        .Eg0_eV = b.Eg0,
        .varshni_alpha_eV_per_K = b.alpha,
        .varshni_beta_K = b.beta,
        .Nc300 = dos_300(b.m_e),
        .Nv300 = dos_300(hole_dos_mass(b.m_hh(), b.m_lh())),
        .electron_affinity_eV = affinity(b.vbo, b.gap_300()),
        .electron_mobility = b.electron,
        .hole_mobility = b.hole,
        .lifetime = {.tau_n0 = b.tau_n, .tau_p0 = b.tau_p, .N_ref = constant_lifetime_N_ref},
        .auger = {.Cn = b.Cn, .Cp = b.Cp},
        .bandgap_narrowing = no_bandgap_narrowing,
        .electron_saturation = {.v_sat_cm_s = b.vsat_n, .beta = 2.0, .temperature_A = b.A_n},
        .hole_saturation = {.v_sat_cm_s = b.vsat_p, .beta = 2.0, .temperature_A = b.A_p},
        .radiative_cm3_s = b.B,
        .ionization = {.donor_eV = 0.0, .acceptor_eV = 0.0, .donor_degeneracy = 2.0, .acceptor_degeneracy = 4.0},
        .richardson = {.electron = 0.0, .hole = 0.0},
    };
    p.band_to_band = {.A_per_cm3_s = 0.0, .B_V_per_cm = 0.0, .direct_gap = b.direct,
                      .electron_mass = b.tunnelling ? b.m_e : 0.0, .hole_mass = b.tunnelling ? b.m_lh() : 0.0};
    p.thermal = {.conductivity_W_cmK = b.k300_W_cmK, .conductivity_exponent = b.k_exponent,
                 .heat_capacity_J_cm3K = b.cp_J_gK * b.rho_g_cm3};
    return p;
}

// One alloy's bowing (amendment 7, Vurgaftman's P(x) = (1 - x) P_A + x P_B - b x (1 - x)).
struct Bowing {
    double gamma_gap, x_gap, m_e, m_hh, m_lh, vbo, C_kappa_W_cmK;
    bool tunnelling;
};

// Binary A (x = 0) to binary B (x = 1), x in (0, 1).
SemiconductorParameters alloy(const Binary& A, const Binary& B, double x, const Bowing& w) {
    const double bow = x * (1.0 - x);
    // Each end's Gamma and X valleys at 300 K, a satellite shifting with the fundamental valley.
    const auto valleys = [](const Binary& b) {
        const double shift = varshni_shift_300(b.alpha, b.beta);
        struct V {
            double gamma, x, m_gamma, m_x;
        };
        const bool x_fundamental = b.satellite_Eg0 > b.Eg0 && !b.direct;  // GaP, AlAs
        if (x_fundamental) return V{b.satellite_Eg0 - shift, b.Eg0 - shift, b.satellite_m_e, b.m_e};
        return V{b.Eg0 - shift, b.satellite_Eg0 > 0.0 ? b.satellite_Eg0 - shift : 0.0, b.m_e, b.satellite_m_e};
    };
    const auto a = valleys(A), c = valleys(B);
    const double gamma = lerp(a.gamma, c.gamma, x) - w.gamma_gap * bow;
    const bool has_x = a.x > 0.0 && c.x > 0.0;
    const double xv = has_x ? lerp(a.x, c.x, x) - w.x_gap * bow : 0.0;
    const bool x_lowest = has_x && xv < gamma;
    const double gap = x_lowest ? xv : gamma;
    const double m_gamma = lerp(a.m_gamma, c.m_gamma, x) - w.m_e * bow;
    const double m_e = x_lowest ? lerp(a.m_x, c.m_x, x) : m_gamma;
    const double m_hh = lerp(A.m_hh(), B.m_hh(), x) - w.m_hh * bow;
    const double m_lh = lerp(A.m_lh(), B.m_lh(), x) - w.m_lh * bow;
    const double vbo = lerp(A.vbo, B.vbo, x) - w.vbo * bow;

    SemiconductorParameters p = make(A);
    p.eps_r = lerp(A.eps, B.eps, x);
    // The set's Varshni is A's; Eg0 puts the 300 K gap at the bowed value, and alloy_gap adds
    // x (d_B - d_A), so the gap's change from 300 K is linear in x.
    p.Eg0_eV = gap + varshni_shift_300(A.alpha, A.beta);
    p.alloy_gap = {.fraction = x, .alpha_a_eV_per_K = A.alpha, .beta_a_K = A.beta, .alpha_b_eV_per_K = B.alpha,
                   .beta_b_K = B.beta};
    p.Nc300 = dos_300(m_e);
    p.Nv300 = dos_300(hole_dos_mass(m_hh, m_lh));
    p.electron_affinity_eV = affinity(vbo, gap);
    p.alloy_mobility = {.fraction = x, .C_mu_electron = 1.0e6, .C_mu_hole = 1.0e6, .electron_b = B.electron,
                        .hole_b = B.hole};
    p.lifetime.tau_n0 = 0.0;  // not given: SRH refuses the alloy (owner ruling F)
    p.lifetime.tau_p0 = 0.0;
    p.auger = {.Cn = lerp(A.Cn, B.Cn, x), .Cp = lerp(A.Cp, B.Cp, x)};
    p.radiative_cm3_s = lerp(A.B, B.B, x);
    p.electron_saturation = {.v_sat_cm_s = lerp(A.vsat_n, B.vsat_n, x), .beta = 2.0,
                             .temperature_A = lerp(A.A_n, B.A_n, x)};
    p.hole_saturation = {.v_sat_cm_s = lerp(A.vsat_p, B.vsat_p, x), .beta = 2.0,
                         .temperature_A = lerp(A.A_p, B.A_p, x)};
    p.band_to_band = {.A_per_cm3_s = 0.0, .B_V_per_cm = 0.0, .direct_gap = !x_lowest,
                      .electron_mass = w.tunnelling ? m_gamma : 0.0, .hole_mass = w.tunnelling ? m_lh : 0.0};
    const double kappa = lerp(A.k300_W_cmK, B.k300_W_cmK, x) - w.C_kappa_W_cmK * bow;
    p.thermal.conductivity_W_cmK = kappa > 0.0 ? kappa : 0.0;  // not positive: no thermal data
    p.thermal.conductivity_exponent = kappa > 0.0 ? 1.17 : 0.0;
    p.thermal.heat_capacity_J_cm3K = lerp(A.cp_J_gK * A.rho_g_cm3, B.cp_J_gK * B.rho_g_cm3, x);
    return p;
}

std::expected<SemiconductorParameters, base::Error> alloy(const char* name, const Binary& A, const Binary& B, double x,
                                                          const Bowing& w) {
    if (!(std::isfinite(x) && x >= 0.0 && x <= 1.0)) {
        return std::unexpected(base::Error{base::ErrorCode::invalid_input,
                                           std::format("the {} In fraction {} lies outside [0, 1]", name, x)});
    }
    if (x == 0.0) return make(A);
    if (x == 1.0) return make(B);
    return alloy(A, B, x, w);
}

}  // namespace

SemiconductorParameters gaas_2001_parameters() { return make(gaas_2001_set); }
SemiconductorParameters inp_parameters() { return make(inp_set); }
SemiconductorParameters inas_parameters() { return make(inas_set); }
SemiconductorParameters gap_parameters() { return make(gap_set); }
SemiconductorParameters alas_parameters() { return make(alas_set); }
SemiconductorParameters gasb_parameters() { return make(gasb_set); }
SemiconductorParameters insb_parameters() { return make(insb_set); }

std::expected<SemiconductorParameters, base::Error> ingaas_parameters(double x) {
    return alloy("InGaAs", gaas_2001_set, inas_set, x, {0.477, 1.40, 0.0091, -0.145, 0.0202, -0.38, 1.41, true});
}

std::expected<SemiconductorParameters, base::Error> inalas_parameters(double x) {
    return alloy("InAlAs", alas_set, inas_set, x, {0.70, 0.0, 0.049, 0.0, 0.0, -0.64, 2.27, false});
}

}  // namespace NiTCAD::physics
