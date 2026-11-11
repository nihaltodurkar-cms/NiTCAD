#include "NiTCAD/verify/reference/trap_admittance.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <utility>

#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/first_integral.hpp"

namespace NiTCAD::verify::reference {

namespace {

using Complex = std::complex<double>;

bool positive_finite(double v) { return std::isfinite(v) && v > 0.0; }

// f_0 = 1 / (1 + e^(tau - u)), from whichever exponential is below 1.
DD fermi(const DD& tau, const DD& u) {
    const DD x = tau - u;
    if (x.hi > 0.0) {
        const DD e = exp(-x);
        return e / (DD(1.0) + e);
    }
    return DD(1.0) / (DD(1.0) + exp(x));
}

// Q_s [C/cm^2] at the surface level u.
DD semiconductor_charge(const TrapAdmittanceInput& in, const NeutralBulk& b, const DD& u) {
    const DD t = u - b.u_b;
    if (t.hi == 0.0 && t.lo == 0.0) return DD(0.0);
    const DD mag = sqrt(DD(2.0) * DD(in.q_C) * DD(in.eps_F_per_cm) * DD(in.thermal_voltage_V) * first_integral_phi(b, u));
    return t.hi > 0.0 || (t.hi == 0.0 && t.lo > 0.0) ? -mag : mag;
}

DD trap_charge(const TrapAdmittanceInput& in, const DD& tau, const DD& u) {
    const DD f = fermi(tau, u);
    const DD qN = DD(in.q_C) * DD(in.density_cm2);
    return in.donor ? qN * (DD(1.0) - f) : -(qN * f);
}

// The neutral substrate's impedance [Ohm cm^2].
Complex bulk_impedance(const TrapOperatingPoint& op, double omega, double neutral_length_cm) {
    return op.sheet_resistance_Ohm_cm * neutral_length_cm / Complex(1.0, omega * op.relaxation_time_s);
}

}  // namespace

std::expected<TrapOperatingPoint, std::string> trap_operating_point(const TrapAdmittanceInput& in) {
    if (!positive_finite(in.q_C) || !positive_finite(in.thermal_voltage_V) || !positive_finite(in.eps_F_per_cm) ||
        !positive_finite(in.ni_cm3) || !std::isfinite(in.net_doping_cm3) || in.net_doping_cm3 == 0.0 ||
        !positive_finite(in.oxide_capacitance_F_per_cm2) || !std::isfinite(in.flat_band_V) ||
        !std::isfinite(in.gate_V) || !std::isfinite(in.density_cm2) || in.density_cm2 < 0.0 ||
        !std::isfinite(in.energy_eV) || !positive_finite(in.capture_n_cm3_s) ||
        !positive_finite(in.capture_p_cm3_s) || !positive_finite(in.mobility_cm2_Vs)) {
        return std::unexpected("trap_admittance: an input is not finite, or not positive where it must be");
    }
    const DD VT(in.thermal_voltage_V), q(in.q_C), eps(in.eps_F_per_cm), cox(in.oxide_capacitance_F_per_cm2);
    const DD tau = DD(in.energy_eV) / VT;
    const double n1 = (DD(in.ni_cm3) * exp(tau)).to_double(), p1 = (DD(in.ni_cm3) * exp(-tau)).to_double();
    if (!positive_finite(n1) || !positive_finite(p1)) {
        return std::unexpected("trap_admittance: the level's n_1 or p_1 is outside the normal range");
    }
    const NeutralBulk b = neutral_bulk(in.ni_cm3, in.net_doping_cm3);
    const DD drive = (DD(in.gate_V) - DD(in.flat_band_V)) / VT;
    // g(u) = (u - u_b) - (Q_s + Q_it) / (C_ox V_T) - (V_G - V_FB) / V_T, rising in u. |u - u_b| <= |drive| + the
    // level's largest shift q N_t / (C_ox V_T) + 1.
    const auto g = [&](const DD& u) {
        return (u - b.u_b) - (semiconductor_charge(in, b, u) + trap_charge(in, tau, u)) / (cox * VT) - drive;
    };
    const double reach = std::abs(drive.hi) + in.q_C * in.density_cm2 / (in.oxide_capacitance_F_per_cm2 * in.thermal_voltage_V) + 1.0;
    const BrentResult root = brent_zero(g, b.u_b - DD(reach), b.u_b + DD(reach), u_dd * (std::abs(b.u_b.hi) + reach));
    if (root.status != BrentStatus::ok) return std::unexpected("trap_admittance: Brent did not converge for u_s");
    TrapOperatingPoint op;
    op.u_s = root.root;
    op.u_b = b.u_b;
    op.brent_evaluations = root.evaluations;
    const DD t = op.u_s - b.u_b;
    op.surface_potential_V = (VT * t).to_double();
    const DD f0 = fermi(tau, op.u_s);
    op.occupancy = f0.to_double();
    op.semiconductor_charge_C_per_cm2 = semiconductor_charge(in, b, op.u_s).to_double();
    op.trap_charge_C_per_cm2 = trap_charge(in, tau, op.u_s).to_double();
    // C_s = sqrt(q eps / (2 V_T)) |Phi'| / sqrt(Phi), Phi' = n_b expm1(t) - p_b expm1(-t); sqrt(q eps (n_b + p_b) /
    // V_T) at t = 0.
    const DD phi = first_integral_phi(b, op.u_s);
    const DD minority_part = in.net_doping_cm3 > 0.0 ? b.p_b * poisson_boltzmann_G(-t) : b.n_b * poisson_boltzmann_G(t);
    if (phi.hi > 0.0) {
        const DD dphi = b.n_b * expm1(t) - b.p_b * expm1(-t);
        op.semiconductor_capacitance_F_per_cm2 = (sqrt(q * eps / (DD(2.0) * VT)) * abs(dphi) / sqrt(phi)).to_double();
        op.p3_minority_charge = (minority_part / phi).to_double();
    } else {
        op.semiconductor_capacitance_F_per_cm2 = sqrt(q * eps * (b.n_b + b.p_b) / VT).to_double();
    }
    op.trap_capacitance_F_per_cm2 = (q * DD(in.density_cm2) * f0 * (DD(1.0) - f0) / VT).to_double();
    const double n0 = (DD(in.ni_cm3) * exp(op.u_s)).to_double(), p0 = (DD(in.ni_cm3) * exp(-op.u_s)).to_double();
    const bool n_type = in.net_doping_cm3 > 0.0;
    op.majority_surface_cm3 = n_type ? n0 : p0;
    op.majority_level_cm3 = n_type ? n1 : p1;
    op.minority_surface_cm3 = n_type ? p0 : n0;
    op.minority_level_cm3 = n_type ? p1 : n1;
    const double c_maj = n_type ? in.capture_n_cm3_s : in.capture_p_cm3_s;
    const double c_min = n_type ? in.capture_p_cm3_s : in.capture_n_cm3_s;
    const double rate = c_maj * (op.majority_surface_cm3 + op.majority_level_cm3);
    op.time_constant_s = 1.0 / rate;
    op.p1_minority_exchange = c_min * (op.minority_surface_cm3 + op.minority_level_cm3) / rate;
    const double N = std::abs(in.net_doping_cm3);
    op.surface_field_V_per_cm = std::abs(op.semiconductor_charge_C_per_cm2) / in.eps_F_per_cm;
    op.depletion_width_cm = std::abs(op.semiconductor_charge_C_per_cm2) / (in.q_C * N);
    op.p2_fermi_drop = op.surface_field_V_per_cm > 0.0
                           ? in.density_cm2 * c_maj / (in.mobility_cm2_Vs * op.surface_field_V_per_cm)
                           : INFINITY;
    op.sheet_resistance_Ohm_cm = 1.0 / (in.q_C * in.mobility_cm2_Vs * N);
    op.relaxation_time_s = in.eps_F_per_cm * op.sheet_resistance_Ohm_cm;
    return op;
}

std::complex<double> trap_admittance(const TrapAdmittanceInput& in, const TrapOperatingPoint& op, double omega,
                                     double neutral_length_cm) {
    const Complex iw(0.0, omega);
    const Complex y_it = iw * op.trap_capacitance_F_per_cm2 / (1.0 + iw * op.time_constant_s);
    const Complex z = 1.0 / (iw * in.oxide_capacitance_F_per_cm2) +
                      1.0 / (iw * op.semiconductor_capacitance_F_per_cm2 + y_it) +
                      bulk_impedance(op, omega, neutral_length_cm);
    return 1.0 / z;
}

std::expected<std::complex<double>, std::string> trap_admittance_direct(const TrapAdmittanceInput& in,
                                                                        const TrapOperatingPoint& op, double omega,
                                                                        double neutral_length_cm) {
    // Unknowns x = (d psi_s, d f, d psi_ox) at a unit gate signal; J = i w C_ox d psi_ox is the gate current.
    //   (1) d psi_ox + d psi_s + Z_b i w C_ox d psi_ox = 1
    //   (2) C_ox d psi_ox - C_s d psi_s - q N_t d f = 0: d Q_G = -d Q_s - d Q_it, d Q_s = -C_s d psi_s, and d Q_it =
    //       -q N_t d f for either trap type (f is the electron occupancy)
    //   (3) (i w + c (m_0 + m_1)) d f - k d psi_s / V_T = 0, the rate equation's majority terms: on n-type
    //       c_n (1 - f_0) d n with d n = n_0 d psi_s / V_T, so k = c_n n_0 (1 - f_0); on p-type -c_p f_0 d p with
    //       d p = -p_0 d psi_s / V_T, so k = c_p p_0 f_0.
    const bool n_type = in.net_doping_cm3 > 0.0;
    const double c = n_type ? in.capture_n_cm3_s : in.capture_p_cm3_s;
    const double f0 = op.occupancy;
    const double k = c * op.majority_surface_cm3 * (n_type ? 1.0 - f0 : f0);
    const Complex iw(0.0, omega);
    const Complex zb = bulk_impedance(op, omega, neutral_length_cm);
    const double cox = in.oxide_capacitance_F_per_cm2;
    std::array<std::array<Complex, 4>, 3> a{{
        {Complex(1.0), Complex(0.0), 1.0 + zb * iw * cox, Complex(1.0)},
        {Complex(-op.semiconductor_capacitance_F_per_cm2), Complex(-in.q_C * in.density_cm2), Complex(cox), Complex(0.0)},
        {Complex(-k / in.thermal_voltage_V), iw + c * (op.majority_surface_cm3 + op.majority_level_cm3), Complex(0.0),
         Complex(0.0)},
    }};
    for (std::size_t col = 0; col < 3; ++col) {
        std::size_t p = col;
        for (std::size_t r = col + 1; r < 3; ++r) {
            if (std::abs(a[r][col]) > std::abs(a[p][col])) p = r;
        }
        if (std::abs(a[p][col]) == 0.0) return std::unexpected("trap_admittance_direct: singular system");
        std::swap(a[p], a[col]);
        for (std::size_t r = 0; r < 3; ++r) {
            if (r == col) continue;
            const Complex m = a[r][col] / a[col][col];
            for (std::size_t j = col; j < 4; ++j) a[r][j] -= m * a[col][j];
        }
    }
    const Complex d_psi_ox = a[2][3] / a[2][2];
    return iw * cox * d_psi_ox;
}

}  // namespace NiTCAD::verify::reference
