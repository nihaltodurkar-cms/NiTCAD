#include "NiTCAD/verify/reference/mos_capacitor.hpp"

#include <cmath>
#include <complex>

#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/complex_step.hpp"

namespace NiTCAD::verify::reference {

namespace {

bool positive_finite(double v) { return std::isfinite(v) && v > 0.0; }

}  // namespace

std::expected<DD, std::string> flat_band_voltage(const FlatBandInput& in) {
    if (!positive_finite(in.work_function_eV) || !std::isfinite(in.electron_affinity_eV) ||
        !positive_finite(in.band_gap_eV) || !positive_finite(in.conduction_dos_cm3) ||
        !positive_finite(in.valence_dos_cm3) || !positive_finite(in.thermal_voltage_V) || !positive_finite(in.ni_cm3) ||
        !std::isfinite(in.net_doping_cm3) || in.net_doping_cm3 == 0.0 || !std::isfinite(in.fixed_charge_cm2) ||
        !positive_finite(in.q_C) || !positive_finite(in.oxide_capacitance_F_per_cm2)) {
        return std::unexpected("flat_band_voltage: an input is not finite, or not positive where it must be");
    }
    const DD VT(in.thermal_voltage_V);
    const DD ec_ei = ldexp(DD(in.band_gap_eV), -1) +
                     ldexp(VT, -1) * log(DD(in.conduction_dos_cm3) / DD(in.valence_dos_cm3));
    const NeutralBulk b = neutral_bulk(in.ni_cm3, in.net_doping_cm3);
    return DD(in.work_function_eV) - DD(in.electron_affinity_eV) - ec_ei + VT * b.u_b -
           DD(in.q_C) * DD(in.fixed_charge_cm2) / DD(in.oxide_capacitance_F_per_cm2);
}

NeutralBulk mos_bulk(const MosInput& in) { return neutral_bulk(in.ni_cm3, in.net_doping_cm3); }

DD mos_debye_cm(const MosInput& in) {
    return sqrt(DD(in.eps_F_per_cm) * DD(in.thermal_voltage_V) / (DD(in.q_C) * DD(std::abs(in.net_doping_cm3))));
}

std::expected<MosPoint, std::string> mos_point(const MosInput& in, double gate_V) {
    if (!positive_finite(in.q_C) || !positive_finite(in.eps_F_per_cm) || !positive_finite(in.thermal_voltage_V) ||
        !positive_finite(in.ni_cm3) || !std::isfinite(in.net_doping_cm3) || in.net_doping_cm3 == 0.0 ||
        !positive_finite(in.oxide_capacitance_F_per_cm2) || !std::isfinite(in.flat_band_V) ||
        !std::isfinite(gate_V)) {
        return std::unexpected("mos_point: an input is not finite, or not positive where it must be");
    }
    const NeutralBulk b = mos_bulk(in);
    const DD VT(in.thermal_voltage_V), Cox(in.oxide_capacitance_F_per_cm2), VFB(in.flat_band_V), VG(gate_V);
    const DD k = sqrt(DD(2.0) * DD(in.q_C) * DD(in.eps_F_per_cm) * VT);  // sqrt(2 q eps V_T)
    const auto charge = [&](const DD& u) {
        const DD t = u - b.u_b;
        const DD mag = k * sqrt(first_integral_phi(b, u));
        return t > DD(0.0) ? -mag : (t < DD(0.0) ? mag : DD(0.0));
    };
    const auto mismatch = [&](const DD& u) { return VFB + VT * (u - b.u_b) - charge(u) / Cox - VG; };
    const DD a = DD(std::abs(gate_V - in.flat_band_V) / in.thermal_voltage_V + 1.0);
    const BrentResult r = brent_zero(mismatch, b.u_b - a, b.u_b + a, u_dd * (std::abs(b.u_b.hi) + a.hi));
    if (r.status != BrentStatus::ok) return std::unexpected("mos_point: Brent did not converge for u_s");
    MosPoint p;
    p.gate_V = VG;
    p.u_s = r.root;
    p.brent_evaluations = r.evaluations;
    const DD t = p.u_s - b.u_b;
    p.band_bending_V = VT * t;
    p.charge_C_per_cm2 = charge(p.u_s);
    if (t.hi == 0.0) {
        p.surface_capacitance = sqrt(DD(in.q_C) * DD(in.eps_F_per_cm) * (b.n_b + b.p_b) / VT);
    } else {
        const DD dphi = b.n_b * expm1(t) - b.p_b * expm1(-t);
        const DD sgn = t > DD(0.0) ? DD(1.0) : DD(-1.0);
        // C_s = -dQ_s / dpsi_s = sgn(t) k Phi' / (2 V_T sqrt Phi).
        p.surface_capacitance = sgn * k * dphi / (DD(2.0) * VT * sqrt(first_integral_phi(b, p.u_s)));
    }
    p.capacitance = p.surface_capacitance * Cox / (p.surface_capacitance + Cox);
    // The second method: dQ_s / du by complex step on Q_s in complex double.
    if (std::abs(t.hi) >= 0.5) {
        const double nb = b.n_b.to_double(), pb = b.p_b.to_double(), ub = b.u_b.to_double(), kk = k.to_double();
        const double sign = t.hi > 0.0 ? -1.0 : 1.0;
        const auto Q = [&](complex_d u) {
            const complex_d tt = u - ub;
            const complex_d phi = nb * (std::exp(tt) - 1.0 - tt) + pb * (std::exp(-tt) - 1.0 + tt);
            return sign * kk * std::sqrt(phi);
        };
        const ComplexStepResult cs = complex_step(Q, p.u_s.to_double());
        if (cs.status == ComplexStepStatus::ok) {
            const double Cs = -cs.derivative / in.thermal_voltage_V;
            p.capacitance_complex_step = Cs * in.oxide_capacitance_F_per_cm2 / (Cs + in.oxide_capacitance_F_per_cm2);
            p.complex_step_evaluated = true;
        }
    }
    return p;
}

}  // namespace NiTCAD::verify::reference
