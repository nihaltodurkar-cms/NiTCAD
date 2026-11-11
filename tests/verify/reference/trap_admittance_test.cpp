// L0 tests of the trap-admittance reference (verify/reference/trap_admittance.hpp; DECISIONS Unit 55 17.3, 17.13
// (a), TVB-22).
//   - without a level (N_t = 0) the operating point is the MOS capacitor reference's (mos_capacitor.hpp): u_s, Q_s
//     and C_s to 1e-13, and Y / (i w) is the series C_ox, C_s;
//   - the gate equation holds at the root, recomputed here in plain double, and a donor-like and an acceptor-like
//     level of one density at one surface level differ in gate bias by q N_t / C_ox exactly;
//   - the level's numbers against their definitions in plain double: f_0, C_T, tau = f_0 / (c_n n_0) on n-type and
//     (1 - f_0) / (c_p p_0) on p-type (not Goetzberger's 1 / (c n_s0): the ratio is f_0), the three predicates;
//   - the circuit against the direct 3 x 3 solve, both substrate types and both trap types, over w tau = 1e-3 ...
//     1e3, to 1e-12 relative;
//   - the branch behind the oxide and the substrate: G_p / w = C_T w tau / (1 + w^2 tau^2), its maximum C_T / 2 at
//     w tau = 1, and the capacitance step C_T from low to high frequency;
//   - the limits: w -> 0 gives C_ox (C_s + C_T) / (C_ox + C_s + C_T), and the substrate adds a conductance w^2 R_b
//     C^2 at low frequency;
//   - bad inputs are refused.
#include "NiTCAD/verify/reference/trap_admittance.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <complex>
#include <limits>

#include "NiTCAD/verify/reference/mos_capacitor.hpp"

namespace ref = NiTCAD::verify::reference;
using Complex = std::complex<double>;

namespace {

constexpr double q = 1.602176634e-19, VT = 0.025851999786435535, eps0 = 8.8541878128e-14;
constexpr double eps_si = 11.7 * eps0, ni = 9.673e9, cox = 3.9 * eps0 / 1e-6;

// An n-type substrate in depletion with an acceptor-like level 0.25 eV above mid-gap.
ref::TrapAdmittanceInput n_type() {
    return {.q_C = q,
            .thermal_voltage_V = VT,
            .eps_F_per_cm = eps_si,
            .ni_cm3 = ni,
            .net_doping_cm3 = 1e16,
            .oxide_capacitance_F_per_cm2 = cox,
            .flat_band_V = -0.2,
            .gate_V = -0.55,
            .donor = false,
            .density_cm2 = 1e10,
            .energy_eV = 0.25,
            .capture_n_cm3_s = 1e-8,
            .capture_p_cm3_s = 4e-9,
            .mobility_cm2_Vs = 1000.0};
}

ref::TrapAdmittanceInput p_type() {
    ref::TrapAdmittanceInput in = n_type();
    in.net_doping_cm3 = -1e16;
    in.flat_band_V = -0.9;
    in.gate_V = -0.55;
    in.donor = true;
    in.energy_eV = -0.25;
    in.mobility_cm2_Vs = 400.0;
    return in;
}

double rel(Complex a, Complex b) { return std::abs(a - b) / std::abs(b); }

}  // namespace

TEST_CASE("trap admittance: without a level the operating point is the MOS capacitor's") {
    ref::TrapAdmittanceInput in = n_type();
    in.density_cm2 = 0.0;
    const auto op = ref::trap_operating_point(in);
    REQUIRE(op.has_value());
    const auto mos = ref::mos_point({.q_C = q,
                                     .eps_F_per_cm = eps_si,
                                     .thermal_voltage_V = VT,
                                     .ni_cm3 = ni,
                                     .net_doping_cm3 = 1e16,
                                     .oxide_capacitance_F_per_cm2 = cox,
                                     .flat_band_V = in.flat_band_V},
                                    in.gate_V);
    REQUIRE(mos.has_value());
    CHECK(std::abs((op->u_s - mos->u_s).to_double()) <= 1e-13 * std::abs(mos->u_s.to_double()));
    CHECK(std::abs(op->semiconductor_charge_C_per_cm2 / mos->charge_C_per_cm2.to_double() - 1.0) <= 1e-13);
    CHECK(std::abs(op->semiconductor_capacitance_F_per_cm2 / mos->surface_capacitance.to_double() - 1.0) <= 1e-13);
    CHECK(op->trap_capacitance_F_per_cm2 == 0.0);
    const double w = 1e3, cs = op->semiconductor_capacitance_F_per_cm2;
    const Complex y = ref::trap_admittance(in, *op, w, 0.0);
    CHECK(rel(y, Complex(0.0, w * cox * cs / (cox + cs))) <= 1e-14);
}

TEST_CASE("trap admittance: the gate equation and the two trap types") {
    for (const auto& in : {n_type(), p_type()}) {
        const auto op = ref::trap_operating_point(in);
        REQUIRE(op.has_value());
        const double residual = in.flat_band_V + op->surface_potential_V -
                                (op->semiconductor_charge_C_per_cm2 + op->trap_charge_C_per_cm2) / cox - in.gate_V;
        INFO("psi_s " << op->surface_potential_V << ", Q_s " << op->semiconductor_charge_C_per_cm2 << ", Q_it "
                      << op->trap_charge_C_per_cm2 << ", f_0 " << op->occupancy);
        CHECK(std::abs(residual) <= 1e-14);
        // Depletion: the surface potential has the sign that repels the majority carriers, Q_s the opposite sign.
        const double sign = in.net_doping_cm3 > 0.0 ? -1.0 : 1.0;
        CHECK(sign * op->surface_potential_V > 0.0);
        CHECK(sign * op->semiconductor_charge_C_per_cm2 < 0.0);
        // The other trap type holds the same surface level at a gate bias shifted by q N_t / C_ox.
        ref::TrapAdmittanceInput other = in;
        other.donor = !in.donor;
        other.gate_V = in.gate_V + (in.donor ? 1.0 : -1.0) * q * in.density_cm2 / cox;
        const auto op2 = ref::trap_operating_point(other);
        REQUIRE(op2.has_value());
        CHECK(std::abs((op2->u_s - op->u_s).to_double()) <= 1e-12);
        const double step = (op2->trap_charge_C_per_cm2 - op->trap_charge_C_per_cm2) / (q * in.density_cm2);
        CHECK(std::abs(step - (in.donor ? -1.0 : 1.0)) <= 1e-12);  // donor-like is the more positive
    }
}

TEST_CASE("trap admittance: the level's numbers against their definitions") {
    for (const auto& in : {n_type(), p_type()}) {
        const auto op = ref::trap_operating_point(in);
        REQUIRE(op.has_value());
        const bool n = in.net_doping_cm3 > 0.0;
        const double us = op->u_s.to_double();
        const double n0 = ni * std::exp(us), p0 = ni * std::exp(-us);
        const double n1 = ni * std::exp(in.energy_eV / VT), p1 = ni * std::exp(-in.energy_eV / VT);
        const double f0 = n0 / (n0 + n1);
        CHECK(std::abs(op->occupancy / f0 - 1.0) <= 1e-13);
        CHECK(std::abs(op->occupancy - p1 / (p0 + p1)) <= 1e-13);  // detailed balance for the holes
        CHECK(std::abs(op->trap_capacitance_F_per_cm2 / (q * in.density_cm2 * f0 * (1.0 - f0) / VT) - 1.0) <= 1e-12);
        const double tau = n ? f0 / (in.capture_n_cm3_s * n0) : (1.0 - f0) / (in.capture_p_cm3_s * p0);
        INFO("f_0 " << f0 << ", tau " << op->time_constant_s << ", P1 " << op->p1_minority_exchange << ", P2 "
                    << op->p2_fermi_drop << ", P3 " << op->p3_minority_charge);
        CHECK(std::abs(op->time_constant_s / tau - 1.0) <= 1e-12);
        // Goetzberger's eq (2), 1 / (c n_s0), is this tau divided by the occupancy factor.
        const double goetzberger = 1.0 / ((n ? in.capture_n_cm3_s * n0 : in.capture_p_cm3_s * p0));
        CHECK(std::abs(op->time_constant_s / goetzberger - (n ? f0 : 1.0 - f0)) <= 1e-12);
        const double p1_pred = n ? in.capture_p_cm3_s * (p0 + p1) / (in.capture_n_cm3_s * (n0 + n1))
                                 : in.capture_n_cm3_s * (n0 + n1) / (in.capture_p_cm3_s * (p0 + p1));
        CHECK(std::abs(op->p1_minority_exchange / p1_pred - 1.0) <= 1e-11);
        const double field = std::abs(op->semiconductor_charge_C_per_cm2) / eps_si;
        CHECK(std::abs(op->surface_field_V_per_cm / field - 1.0) <= 1e-14);
        const double c = n ? in.capture_n_cm3_s : in.capture_p_cm3_s;
        CHECK(std::abs(op->p2_fermi_drop / (in.density_cm2 * c / (in.mobility_cm2_Vs * field)) - 1.0) <= 1e-13);
        CHECK(op->p3_minority_charge > 0.0);
        CHECK(op->p3_minority_charge < 1e-6);  // depletion: the minority carriers hold almost none of the charge
        CHECK(std::abs(op->depletion_width_cm * q * 1e16 / std::abs(op->semiconductor_charge_C_per_cm2) - 1.0) <= 1e-14);
        CHECK(std::abs(op->relaxation_time_s * q * in.mobility_cm2_Vs * 1e16 / eps_si - 1.0) <= 1e-14);
    }
}

TEST_CASE("trap admittance: the circuit agrees with the direct solve") {
    for (ref::TrapAdmittanceInput in : {n_type(), p_type()}) {
        for (const bool donor : {false, true}) {
            in.donor = donor;
            const auto op = ref::trap_operating_point(in);
            REQUIRE(op.has_value());
            for (int k = -6; k <= 6; ++k) {
                const double w = std::pow(10.0, 0.5 * k) / op->time_constant_s;
                for (const double Ln : {0.0, 2e-4}) {
                    const Complex a = ref::trap_admittance(in, *op, w, Ln);
                    const auto b = ref::trap_admittance_direct(in, *op, w, Ln);
                    REQUIRE(b.has_value());
                    INFO("w tau " << w * op->time_constant_s << ", L_n " << Ln << ": " << a.real() << " + i "
                                  << a.imag() << " against " << b->real() << " + i " << b->imag());
                    CHECK(rel(a, *b) <= 1e-12);
                    CHECK(a.real() >= 0.0);  // passive
                }
            }
        }
    }
}

TEST_CASE("trap admittance: the conductance peak and the capacitance step behind the oxide") {
    const ref::TrapAdmittanceInput in = n_type();
    const auto op = ref::trap_operating_point(in);
    REQUIRE(op.has_value());
    const double tau = op->time_constant_s, CT = op->trap_capacitance_F_per_cm2, cs = op->semiconductor_capacitance_F_per_cm2;
    const auto behind = [&](double w) {  // with no substrate resistance
        const Complex y = ref::trap_admittance(in, *op, w, 0.0);
        return 1.0 / (1.0 / y - 1.0 / Complex(0.0, w * cox));
    };
    double peak = 0.0, at = 0.0;
    for (int k = -200; k <= 200; ++k) {
        const double wt = std::pow(10.0, 0.01 * k), w = wt / tau;
        const Complex ys = behind(w);
        const double g = ys.real() / w;
        CHECK(std::abs(g - CT * wt / (1.0 + wt * wt)) <= 1e-11 * CT);
        CHECK(std::abs(ys.imag() / w - (cs + CT / (1.0 + wt * wt))) <= 1e-11 * (cs + CT));
        if (g > peak) {
            peak = g;
            at = wt;
        }
    }
    CHECK(at == 1.0);
    CHECK(std::abs(peak / (0.5 * CT) - 1.0) <= 1e-11);
    const double low = (behind(1e-6 / tau)).imag() / (1e-6 / tau), high = (behind(1e6 / tau)).imag() / (1e6 / tau);
    CHECK(std::abs((low - high) / CT - 1.0) <= 1e-9);
}

TEST_CASE("trap admittance: the low-frequency limit and the substrate's loss") {
    const ref::TrapAdmittanceInput in = n_type();
    const auto op = ref::trap_operating_point(in);
    REQUIRE(op.has_value());
    const double tau = op->time_constant_s, CT = op->trap_capacitance_F_per_cm2, cs = op->semiconductor_capacitance_F_per_cm2;
    const double w = 1e-6 / tau;
    const double C_lf = cox * (cs + CT) / (cox + cs + CT);
    CHECK(std::abs(ref::trap_admittance(in, *op, w, 0.0).imag() / (w * C_lf) - 1.0) <= 1e-10);
    // With a neutral length the series resistance R_b adds Re Y = w^2 R_b C^2 (far below the relaxation rate).
    const double Ln = 1e-3, Rb = op->sheet_resistance_Ohm_cm * Ln;
    const double extra = ref::trap_admittance(in, *op, w, Ln).real() - ref::trap_admittance(in, *op, w, 0.0).real();
    INFO("R_b " << Rb << " Ohm cm^2, extra Re Y " << extra << ", w^2 R_b C^2 " << w * w * Rb * C_lf * C_lf);
    CHECK(std::abs(extra / (w * w * Rb * C_lf * C_lf) - 1.0) <= 1e-6);
    // At w tau_d = 1 a long substrate dominates the impedance: Y = (1 + i w tau_d) / R_b = (1 + i) / R_b.
    const double wd = 1.0 / op->relaxation_time_s, Rlong = op->sheet_resistance_Ohm_cm * 1e4;
    const Complex y = ref::trap_admittance(in, *op, wd, 1e4);
    CHECK(std::abs(y.real() * Rlong - 1.0) <= 1e-6);
    CHECK(std::abs(y.imag() * Rlong - 1.0) <= 1e-6);
}

TEST_CASE("trap admittance: bad inputs are refused") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const auto bad = [](const ref::TrapAdmittanceInput& in) { return !ref::trap_operating_point(in).has_value(); };
    CHECK_FALSE(bad(n_type()));
    for (double ref::TrapAdmittanceInput::* m :
         {&ref::TrapAdmittanceInput::q_C, &ref::TrapAdmittanceInput::thermal_voltage_V,
          &ref::TrapAdmittanceInput::eps_F_per_cm, &ref::TrapAdmittanceInput::ni_cm3,
          &ref::TrapAdmittanceInput::net_doping_cm3, &ref::TrapAdmittanceInput::oxide_capacitance_F_per_cm2,
          &ref::TrapAdmittanceInput::capture_n_cm3_s, &ref::TrapAdmittanceInput::capture_p_cm3_s,
          &ref::TrapAdmittanceInput::mobility_cm2_Vs}) {
        ref::TrapAdmittanceInput in = n_type();
        in.*m = 0.0;
        CHECK(bad(in));
        in.*m = nan;
        CHECK(bad(in));
    }
    ref::TrapAdmittanceInput in = n_type();
    in.density_cm2 = -1.0;
    CHECK(bad(in));
    in = n_type();
    in.energy_eV = 40.0;  // n_1 overflows
    CHECK(bad(in));
    in = n_type();
    in.gate_V = nan;
    CHECK(bad(in));
}
