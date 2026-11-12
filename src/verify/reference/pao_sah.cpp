#include "NiTCAD/verify/reference/pao_sah.hpp"

#include <algorithm>
#include <cmath>

#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/first_integral.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"

namespace NiTCAD::verify::reference {

namespace {

constexpr double inner_rel = 1e-24, outer_rel = 1e-22;

bool positive_finite(double v) { return std::isfinite(v) && v > 0.0; }

std::expected<void, std::string> check(const PaoSahInput& in, double gate_V, double drain_V) {
    if (!positive_finite(in.q_C) || !positive_finite(in.eps_F_per_cm) || !positive_finite(in.thermal_voltage_V) ||
        !positive_finite(in.ni_cm3) || !positive_finite(in.acceptors_cm3) ||
        !positive_finite(in.oxide_capacitance_F_per_cm2) || !std::isfinite(in.flat_band_V) ||
        !positive_finite(in.mobility_cm2_Vs) || (in.form != 1 && in.form != 2) ||
        !(in.cutoff_fraction > 0.0 && in.cutoff_fraction <= 1.0)) {
        return std::unexpected("pao_sah: an input is not finite, or out of range");
    }
    if (!std::isfinite(gate_V) || !(gate_V > in.flat_band_V)) return std::unexpected("pao_sah: V_G <= V_FB");
    if (!std::isfinite(drain_V) || drain_V < 0.0) return std::unexpected("pao_sah: V_D < 0");
    return {};
}

struct Model {
    const PaoSahInput& in;
    NeutralBulk bulk;
    DD VT, scale;  // E^2 = scale [...]
    explicit Model(const PaoSahInput& i)
        : in(i),
          bulk(neutral_bulk(i.ni_cm3, -i.acceptors_cm3)),
          VT(i.thermal_voltage_V),
          scale(DD(2.0) * DD(i.q_C) * DD(i.thermal_voltage_V) / DD(i.eps_F_per_cm)) {}
    DD psi_B() const { return VT * log(bulk.p_b / DD(in.ni_cm3)); }
    // E(psi, V) [V/cm], 0 where form 1's E^2 is negative.
    DD field(const DD& psi, const DD& V) const {
        const DD u = psi / VT, v = V / VT;
        const DD holes = bulk.p_b * (expm1(-u) + u);
        const DD ev = exp(-v);
        const DD electrons = in.form == 1 ? bulk.n_b * (ev * expm1(u) - u) : bulk.n_b * (ev * (expm1(u) - u));
        const DD e2 = scale * (holes + electrons);
        return e2.hi > 0.0 ? sqrt(e2) : DD(0.0);
    }
    DD n(const DD& psi, const DD& V) const { return bulk.n_b * exp((psi - V) / VT); }
};

std::expected<DD, std::string> surface(const Model& m, double gate_V, const DD& V) {
    const DD drive = DD(gate_V) - DD(m.in.flat_band_V);
    const DD k = DD(m.in.eps_F_per_cm) / DD(m.in.oxide_capacitance_F_per_cm2);
    const auto f = [&](const DD& psi) { return drive - psi - k * m.field(psi, V); };
    const BrentResult b = brent_zero(f, DD(0.0), drive, 1e-30 * drive.hi);
    if (b.status != BrentStatus::ok) return std::unexpected("pao_sah: the surface potential did not converge");
    return b.root;
}

bool quadrature_ok(QuadratureStatus s) { return s == QuadratureStatus::ok || s == QuadratureStatus::roundoff; }

std::expected<InversionCharge, std::string> inversion(const Model& m, double gate_V, const DD& V) {
    const PaoSahInput& in = m.in;
    const auto ps = surface(m, gate_V, V);
    if (!ps) return std::unexpected(ps.error());
    const DD lower = DD(in.cutoff_fraction) * m.psi_B();
    InversionCharge out;
    out.surface_V = *ps;
    if (!(lower < *ps)) {  // the surface below the cutoff: no inversion charge is counted
        out.charge = DD(0.0);
        return out;
    }
    const AdaptiveResult r =
        qag15([&](const DD& psi) { return m.n(psi, V) / m.field(psi, V); }, lower, *ps, 0.0, inner_rel);
    if (!quadrature_ok(r.status)) return std::unexpected("pao_sah: the inner integral failed");
    out.charge = -DD(in.q_C) * r.value;
    out.quadrature_error = r.value.hi != 0.0 ? r.abs_error / std::abs(r.value.hi) : 0.0;
    return out;
}

}  // namespace

std::expected<InversionCharge, std::string> pao_sah_inversion_charge(const PaoSahInput& in, double gate_V,
                                                                     double channel_V) {
    if (auto c = check(in, gate_V, channel_V); !c) return std::unexpected(c.error());
    return inversion(Model(in), gate_V, DD(channel_V));
}

std::expected<PaoSahPoint, std::string> pao_sah(const PaoSahInput& in, double gate_V, double drain_V) {
    if (auto c = check(in, gate_V, drain_V); !c) return std::unexpected(c.error());
    const Model m(in);
    PaoSahPoint out;
    out.bulk_potential_V = m.psi_B();
    const auto s0 = surface(m, gate_V, DD(0.0));
    const auto sL = surface(m, gate_V, DD(drain_V));
    if (!s0 || !sL) return std::unexpected(!s0 ? s0.error() : sL.error());
    out.surface_source_V = *s0;
    out.surface_drain_V = *sL;
    if (drain_V == 0.0) return out;
    double inner_worst = 0.0;
    bool failed = false;
    const auto integrand = [&](const DD& V) {
        const auto q = inversion(m, gate_V, V);
        if (!q) {
            failed = true;
            return DD(0.0);
        }
        inner_worst = std::max(inner_worst, q->quadrature_error);
        return -q->charge;
    };
    const AdaptiveResult r = qag15(integrand, DD(0.0), DD(drain_V), 0.0, outer_rel);
    if (failed || !quadrature_ok(r.status)) return std::unexpected("pao_sah: the outer integral failed");
    out.current = DD(in.mobility_cm2_Vs) * r.value;
    out.quadrature_error = (r.value.hi != 0.0 ? r.abs_error / std::abs(r.value.hi) : 0.0) + inner_worst;
    return out;
}

std::expected<double, std::string> pao_sah_conductance_identity(const PaoSahInput& in, double gate_V,
                                                                double drain_V, double h_V) {
    if (!(h_V > 0.0 && 2.0 * h_V < drain_V)) return std::unexpected("pao_sah: the difference step is out of range");
    // The fourth-order central difference (its truncation about (h / V_T)^4 / 30 relative).
    const auto a = pao_sah(in, gate_V, drain_V - 2.0 * h_V), b = pao_sah(in, gate_V, drain_V - h_V);
    const auto c = pao_sah(in, gate_V, drain_V + h_V), d = pao_sah(in, gate_V, drain_V + 2.0 * h_V);
    const auto q = pao_sah_inversion_charge(in, gate_V, drain_V);
    if (!a || !b || !c || !d || !q) return std::unexpected("pao_sah: a point of the difference failed");
    const DD gd = (a->current - DD(8.0) * b->current + DD(8.0) * c->current - d->current) / (DD(12.0) * DD(h_V));
    const DD expected = DD(in.mobility_cm2_Vs) * -q->charge;
    return std::abs(((gd - expected) / expected).to_double());
}

std::expected<ChargeSheet, std::string> charge_sheet(const PaoSahInput& in, double gate_V, double drain_V) {
    PaoSahInput two = in;
    two.form = 2;
    if (auto c = check(two, gate_V, drain_V); !c) return std::unexpected(c.error());
    const Model m(two);
    const auto s0 = surface(m, gate_V, DD(0.0));
    const auto sL = surface(m, gate_V, DD(drain_V));
    if (!s0 || !sL) return std::unexpected(!s0 ? s0.error() : sL.error());
    const DD cox(in.oxide_capacitance_F_per_cm2), drive = DD(gate_V) - DD(in.flat_band_V);
    const DD gamma = sqrt(DD(2.0) * DD(in.q_C) * DD(in.eps_F_per_cm) * DD(in.acceptors_cm3));
    const DD VT(in.thermal_voltage_V);
    if (!(*s0 > VT) || !(*sL > VT)) return std::unexpected("charge_sheet: psi_s <= V_T");
    const auto QI = [&](const DD& psi) { return -cox * (drive - psi) + gamma * sqrt(psi - VT); };
    // int -Q_I dpsi = cox (drive psi - psi^2 / 2) - gamma (2/3) (psi - V_T)^{3/2}
    const auto F = [&](const DD& psi) {
        const DD t = psi - VT;
        return cox * (drive * psi - psi * psi / DD(2.0)) - gamma * DD(2.0) / DD(3.0) * t * sqrt(t);
    };
    const DD drift = F(*sL) - F(*s0);
    const AdaptiveResult r = qag15([&](const DD& psi) { return -QI(psi); }, *s0, *sL, 0.0, 1e-26);
    if (!quadrature_ok(r.status)) return std::unexpected("charge_sheet: the quadrature failed");
    ChargeSheet out;
    out.current = DD(in.mobility_cm2_Vs) * (drift + VT * (QI(*sL) - QI(*s0)));
    out.second_method = drift.hi != 0.0 ? std::abs(((r.value - drift) / drift).to_double()) : 0.0;
    return out;
}

}  // namespace NiTCAD::verify::reference
