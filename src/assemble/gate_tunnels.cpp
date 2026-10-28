// Gate tunnelling in the drift-diffusion system (Unit 51; DECISIONS.md GD1-GD8): the tunnels'
// geometry from the device, and each tunnel's sinks and partials at a state.
//
// Energies in eV on the band diagram's scale (band_edges.hpp times V_T): the equilibrium Fermi level
// at 0, the vacuum level depth - V_T psi (depth the reference material's intrinsic-level depth), so
//     E_c = V_T (g_n - psi - s),  E_v = -V_T (g_p + psi + s),
//     E_Fn = V_T (ln(n / n_ie) - ln gamma_n - psi - s),  E_Fp = V_T (ln(p / n_ie)... negated),
// an electrode's Fermi level -V (its node's psi_E - electrode_potential, times -V_T) and a lumped
// gate's -V_G. The insulator's conduction edge is the vacuum level less its affinity, its valence
// edge E_g,ox lower (GD1). A hole process runs on the mirrored scale (physics/gate_tunnelling.hpp).
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <variant>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/gate.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/gate_tunnelling.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/thermionic_emission.hpp"
#include "scaled_device.hpp"

namespace NiTCAD::assemble {

namespace {

base::Error refusal(std::string message) {
    return base::Error{base::ErrorCode::invalid_input, std::move(message), {}};
}

// The band data a tunnel through an insulator uses: its own, or SiO2's for the built-in SiO2 set
// (GD5); nullopt without any.
std::optional<physics::InsulatorParameters> band_data(const physics::InsulatorParameters& p) {
    if (physics::has_band_data(p)) return p;
    if (p == physics::silicon_dioxide_parameters) {
        return physics::silicon_dioxide_tunnelling_parameters;
    }
    return std::nullopt;
}

}  // namespace

std::optional<base::Error> DriftDiffusion::make_gate_tunnels(const device::Device& device,
                                                             const Scaling& scaling) {
    const auto contacts = device.contacts();
    if (std::ranges::none_of(contacts, [](const device::Contact& c) { return c.tunnelling != 0; })) {
        return std::nullopt;
    }
    const mesh::Mesh& mesh = device.mesh();
    const int D = mesh.dimension();
    const double T = scaling.temperature_K;
    const double flow = 1.0 / (base::q_C * scaling.Ns * std::pow(scaling.L_D, D - 2) * scaling.D0);
    const double depth =
        physics::intrinsic_level_depth_eV(device.material(device.reference_node()), T);
    const auto richardson = [T](const physics::Semiconductor& m, physics::Carrier c) {
        const double N = c == physics::Carrier::electron ? physics::conduction_band_dos(m, T)
                                                         : physics::valence_band_dos(m, T);
        return base::q_C * physics::emission_velocity_cm_s(m, c, T) * N / (T * T);
    };
    const auto insulator_of = [&](std::size_t v) -> const physics::InsulatorParameters& {
        const device::Region& r =
            device.regions()[static_cast<std::size_t>(device.node_region()[v])];
        return std::get<physics::Insulator>(r.material).parameters();
    };
    const auto permittivity = [&](std::size_t v) {
        return device.is_insulator(static_cast<mesh::NodeId>(v))
                   ? insulator_of(v).eps_r
                   : device.material(static_cast<mesh::NodeId>(v)).parameters().eps_r;
    };
    const auto piece = [](std::size_t p0, std::size_t p1, double length_cm,
                          const physics::InsulatorParameters& b) {
        return TunnelPiece{p0, p1, length_cm, b.electron_affinity_eV, b.band_gap_eV,
                           b.tunnel_mass_n, b.tunnel_mass_p, physics::evb_mass(b)};
    };
    const auto base_tunnel = [&](std::size_t node, std::size_t c) {
        const physics::Semiconductor& m = device.material(static_cast<mesh::NodeId>(node));
        GateTunnel t{};
        t.node = node;
        t.contact = static_cast<std::int32_t>(c);
        t.processes = contacts[c].tunnelling;
        t.richardson_n = richardson(m, physics::Carrier::electron);
        t.richardson_p = richardson(m, physics::Carrier::hole);
        t.depth_eV = depth;
        return t;
    };

    for (std::size_t c = 0; c < contacts.size(); ++c) {
        const device::Contact& contact = contacts[c];
        if (contact.tunnelling == 0) continue;
        if (contact.kind == device::ContactKind::gate) {
            // Lumped (GD2): one piece of t_ox from the node to the gate, the stack's band data or
            // SiO2's (the device has checked a stack without its own is SiO2's permittivity).
            const device::GateStack& g = contact.gate;
            physics::InsulatorParameters b = physics::silicon_dioxide_tunnelling_parameters;
            if (g.oxide_band_gap_eV > 0.0) {
                b.electron_affinity_eV = g.oxide_affinity_eV;
                b.band_gap_eV = g.oxide_band_gap_eV;
                b.tunnel_mass_n = g.oxide_mass_n;
                b.tunnel_mass_p = g.oxide_mass_p;
                b.tunnel_mass_evb = g.oxide_mass_evb;
            }
            const mesh::BoundaryPatch* patch = mesh.find_boundary(g.boundary);
            NITCAD_EXPECTS(patch != nullptr);  // checked by the device
            for (const mesh::NodeId v : contact.nodes) {
                const auto it = std::ranges::lower_bound(patch->nodes, v);
                const auto k = static_cast<std::size_t>(it - patch->nodes.begin());
                const auto i = static_cast<std::size_t>(v);
                const physics::Semiconductor& m = device.material(v);
                GateTunnel t = base_tunnel(i, c);
                t.lumped = true;
                t.electrode = g.electrode;
                t.work_function_eV = gate_work_function_eV(g, m, T);
                t.gate_gap_eV = physics::band_gap_eV(m, T);
                t.coupling = patch->areas[k] * flow;
                t.points = {{0, 0, 1.0}, {0, 0, 1.0}};  // the second is the gate's own
                t.pieces = {piece(0, 1, g.oxide_thickness_cm, b)};
                gate_tunnels_.push_back(std::move(t));
            }
            continue;
        }
        // An electrode (GD2): straight grid lines from each semiconductor node through the
        // insulator to one of its nodes.
        const mesh::TensorCells* cells = device.cells();
        if (cells == nullptr) {
            return refusal("electrode '" + contact.name +
                           "': gate tunnelling through a meshed insulator needs a tensor grid");
        }
        std::vector<char> mine(mesh.node_count(), 0), other(mesh.node_count(), 0);
        for (std::size_t q = 0; q < contacts.size(); ++q) {
            for (const mesh::NodeId v : contacts[q].nodes) {
                (q == c ? mine : other)[static_cast<std::size_t>(v)] = 1;
            }
        }
        const double phi_m = electrode_work_function_eV(contact.electrode, T);
        const double si_gap = physics::band_gap_eV(physics::silicon(), T);
        for (const mesh::Edge& e : mesh.edges()) {
            const bool first_insulator = device.is_insulator(e.first);
            if (first_insulator == device.is_insulator(e.second)) continue;
            const auto i = static_cast<std::size_t>(first_insulator ? e.second : e.first);
            const auto j = static_cast<std::size_t>(first_insulator ? e.first : e.second);
            if (other[j] != 0 && mine[j] == 0) continue;  // another contact's electrode node
            const mesh::TensorCells::Index ii = cells->index(static_cast<mesh::NodeId>(i));
            const mesh::TensorCells::Index ij = cells->index(static_cast<mesh::NodeId>(j));
            int axis = 0;
            for (int a = 0; a < D; ++a) {
                if (ii[a] != ij[a]) axis = a;
            }
            const bool up = ij[axis] > ii[axis];
            std::vector<std::size_t> path{j};
            bool reached = mine[j] != 0;
            mesh::TensorCells::Index at = ij;
            while (!reached) {
                if (up ? at[axis] + 1 >= cells->size(axis) : at[axis] == 0) break;
                at[axis] = up ? at[axis] + 1 : at[axis] - 1;
                const auto next = static_cast<std::size_t>(cells->node(at));
                if (!device.is_insulator(static_cast<mesh::NodeId>(next)) ||
                    (other[next] != 0 && mine[next] == 0)) {
                    break;
                }
                path.push_back(next);
                reached = mine[next] != 0;
            }
            if (!reached) continue;
            for (const std::size_t v : path) {
                if (!band_data(insulator_of(v))) {
                    return refusal("electrode '" + contact.name +
                                   "': a tunnelling path crosses an insulator without band data "
                                   "(affinity, band gap, tunnelling masses)");
                }
            }
            GateTunnel t = base_tunnel(i, c);
            t.lumped = false;
            t.electrode = contact.electrode.kind;
            t.work_function_eV = phi_m;
            t.gate_gap_eV = si_gap;
            t.coupling = e.coupling_area * flow;
            t.nodes = path;
            // Points and pieces: the interface at the first edge's midpoint (the semiconductor's
            // half is not barrier), then each edge whole within one material or split at its
            // midpoint, where the potential is the permittivity-weighted one of the box method.
            const auto column = [](std::size_t k) { return 3 + k; };
            const auto midpoint = [&](std::size_t ca, std::size_t va, std::size_t cb,
                                      std::size_t vb) {
                const double ea = permittivity(va), eb = permittivity(vb);
                t.points.push_back({ca, cb, ea / (ea + eb)});
                return t.points.size() - 1;
            };
            const auto point_of = [&](std::size_t k) {
                t.points.push_back({column(k), column(k), 1.0});
                return t.points.size() - 1;
            };
            const double first_length =
                std::abs(cells->axis(axis)[ij[axis]] - cells->axis(axis)[ii[axis]]);
            std::size_t from = midpoint(0, i, column(0), j);
            std::size_t to = point_of(0);
            t.pieces.push_back(piece(from, to, 0.5 * first_length, *band_data(insulator_of(j))));
            for (std::size_t k = 0; k + 1 < path.size(); ++k) {
                const std::size_t u = path[k], v = path[k + 1];
                const mesh::TensorCells::Index iu = cells->index(static_cast<mesh::NodeId>(u));
                const mesh::TensorCells::Index iv = cells->index(static_cast<mesh::NodeId>(v));
                const double L = std::abs(cells->axis(axis)[iv[axis]] - cells->axis(axis)[iu[axis]]);
                const physics::InsulatorParameters bu = *band_data(insulator_of(u));
                const physics::InsulatorParameters bv = *band_data(insulator_of(v));
                from = to;
                if (insulator_of(u) == insulator_of(v)) {
                    to = point_of(k + 1);
                    t.pieces.push_back(piece(from, to, L, bu));
                } else {
                    const std::size_t mid = midpoint(column(k), u, column(k + 1), v);
                    t.pieces.push_back(piece(from, mid, 0.5 * L, bu));
                    to = point_of(k + 1);
                    t.pieces.push_back(piece(mid, to, 0.5 * L, bv));
                }
            }
            gate_tunnels_.push_back(std::move(t));
        }
    }
    return std::nullopt;
}

void DriftDiffusion::locate_gate_tunnels() {
    for (GateTunnel& t : gate_tunnels_) {
        t.pos_n.clear();
        t.pos_p.clear();
        for (const std::size_t v : t.nodes) {
            t.pos_n.push_back(detail::position(pattern_, m_ * t.node + 1, m_ * v));
            t.pos_p.push_back(detail::position(pattern_, m_ * t.node + 2, m_ * v));
        }
    }
}

void DriftDiffusion::gate_tunnel_flux(const GateTunnel& t, std::span<const double> x,
                                      std::span<const NodeDegeneracy> g,
                                      GateTunnelFlux& out) const {
    const std::size_t i = t.node;
    const std::size_t K = 3 + t.nodes.size() + (t.lumped ? 1 : 0);
    const std::size_t bias = 3 + t.nodes.size();  // a lumped gate's bias column
    out.T_n = out.T_p = 0.0;
    out.J_ecb = out.J_hvb = out.J_evb = 0.0;
    out.d_n.assign(K, 0.0);
    out.d_p.assign(K, 0.0);

    // A quantity linear in the local columns near x: its value and gradient.
    struct Q {
        double v;
        std::vector<double> d;
    };
    const auto zero = [K](double v) { return Q{v, std::vector<double>(K, 0.0)}; };
    const double V_T = V_T_;
    const double psi = x[m_ * i], n = x[m_ * i + 1], p = x[m_ * i + 2];
    const double s = band_shift_[i], nie = n_ie_[i];
    const double lgn = g.empty() ? 0.0 : g[i].n.log_gamma;
    const double dgn = g.empty() ? 0.0 : g[i].n.d_density;
    const double lgp = g.empty() ? 0.0 : g[i].p.log_gamma;
    const double dgp = g.empty() ? 0.0 : g[i].p.d_density;
    Q Ec = zero(V_T * (log_dos_n_[i] - psi - s));
    Ec.d[0] = -V_T;
    Q Ev = zero(-V_T * (log_dos_p_[i] + psi + s));
    Ev.d[0] = -V_T;
    Q Fn = zero(V_T * (std::log(n / nie) - lgn - psi - s));
    Fn.d[0] = -V_T;
    Fn.d[1] = V_T * (1.0 / n - dgn);
    Q Fp = zero(-V_T * (psi + s + std::log(p / nie) - lgp));
    Fp.d[0] = -V_T;
    Fp.d[2] = -V_T * (1.0 / p - dgp);
    // The electrode's Fermi level.
    Q FG = zero(0.0);
    if (t.lumped) {
        FG.v = -V_T * (gates_.potential(i) + gates_.term(i).offset);  // -V_G
        FG.d[bias] = -1.0;
    } else {
        const std::size_t e = t.nodes.back();
        FG.v = -V_T * (x[m_ * e] - electrode_potential_[e]);
        FG.d[3 + t.nodes.size() - 1] = -V_T;
    }
    const auto column_psi = [&](std::size_t col) {
        return col == 0 ? psi : x[m_ * t.nodes[col - 3]];
    };
    // Vacuum level of each point.
    std::vector<Q> vac;
    vac.reserve(t.points.size());
    for (std::size_t k = 0; k < t.points.size(); ++k) {
        if (t.lumped && k == 1) {
            Q q = zero(FG.v + t.work_function_eV);
            q.d = FG.d;
            vac.push_back(std::move(q));
            continue;
        }
        const TunnelPoint& pt = t.points[k];
        Q q = zero(t.depth_eV - V_T * (pt.w * column_psi(pt.a) + (1.0 - pt.w) * column_psi(pt.b)));
        q.d[pt.a] -= V_T * pt.w;
        q.d[pt.b] -= V_T * (1.0 - pt.w);
        vac.push_back(std::move(q));
    }
    // The electrode's band edges: a metal has none (+-infinity), polysilicon is degenerate.
    constexpr double inf = std::numeric_limits<double>::infinity();
    const bool metal = t.electrode == device::GateElectrode::metal;
    Q EcG = FG, EvG = FG;
    if (metal) {
        EcG.v = -inf;
        EvG.v = inf;
    } else if (t.electrode == device::GateElectrode::n_poly) {
        EvG.v -= t.gate_gap_eV;
    } else {
        EcG.v += t.gate_gap_eV;
    }
    const auto negate = [](Q q) {
        q.v = -q.v;
        for (double& d : q.d) d = -d;
        return q;
    };
    const auto larger = [](const Q& a, const Q& b) -> const Q& { return b.v > a.v ? b : a; };

    enum class Band { conduction, valence };
    std::vector<physics::BarrierSegment> segments(t.pieces.size());
    std::vector<double> dW(2 * t.pieces.size());
    // One process: the barrier on its scale (mirrored for holes), the range [low, high] (high
    // infinite: the barrier top), the Fermi levels; adds J to `value` and its gradient times `sign`
    // and `coupling` to `grad`.
    const auto run = [&](Band band, bool mirrored, const Q& low,
                         const Q* high, const Q& F1, const Q& F2, double richardson, double sign,
                         double& value, std::vector<double>& grad, int mass_kind) {
        for (std::size_t k = 0; k < t.pieces.size(); ++k) {
            const TunnelPiece& pc = t.pieces[k];
            const double offset = band == Band::conduction ? pc.affinity_eV : pc.affinity_eV + pc.gap_eV;
            const double w0 = vac[pc.p0].v - offset, w1 = vac[pc.p1].v - offset;
            const double m = mass_kind == 0 ? pc.mass_n : mass_kind == 1 ? pc.mass_p : pc.mass_evb;
            segments[k] = mirrored ? physics::BarrierSegment{-w0, -w1, pc.length_cm, m}
                                   : physics::BarrierSegment{w0, w1, pc.length_cm, m};
        }
        std::ranges::fill(dW, 0.0);
        const double h = high == nullptr ? inf : high->v;
        const physics::TunnelCurrent J = physics::tsu_esaki_current(
            segments, {low.v, h}, F1.v, F2.v, temperature_K_, richardson, dW);
        if (J.J == 0.0 && J.d_F1 == 0.0 && J.d_F2 == 0.0) return 0.0;
        const double c = sign * t.coupling;
        value += c * J.J;
        const double ws = mirrored ? -1.0 : 1.0;  // d W / d vacuum
        for (std::size_t col = 0; col < K; ++col) {
            double d = J.d_F1 * F1.d[col] + J.d_F2 * F2.d[col] + J.d_E_low * low.d[col];
            if (high != nullptr) d += J.d_E_high * high->d[col];
            for (std::size_t k = 0; k < t.pieces.size(); ++k) {
                d += ws * (dW[2 * k] * vac[t.pieces[k].p0].d[col] +
                           dW[2 * k + 1] * vac[t.pieces[k].p1].d[col]);
            }
            grad[col] += c * d;
        }
        return J.J;
    };

    if (device::tunnels(t.processes, device::TunnelProcess::ecb)) {
        // Electrons from the conduction band to the electrode's (GD4).
        const Q& low = metal ? Ec : larger(Ec, EcG);
        out.J_ecb = run(Band::conduction, false, low, nullptr, Fn, FG, t.richardson_n, 1.0,
                        out.T_n, out.d_n, 0);
    }
    if (device::tunnels(t.processes, device::TunnelProcess::hvb)) {
        // Holes, on the mirrored scale: from -min(E_v, E_v,G) up, levels -E_Fp and -E_FG.
        const Q low = metal ? negate(Ev) : larger(negate(Ev), negate(EvG));
        out.J_hvb = run(Band::valence, true, low, nullptr, negate(Fp), negate(FG),
                        t.richardson_p, 1.0, out.T_p, out.d_p, 1);
    }
    if (device::tunnels(t.processes, device::TunnelProcess::evb) && !metal) {
        // Valence electrons of the semiconductor into the electrode's conduction band (a hole
        // left behind: T_p falls), and the electrode's valence electrons into the
        // semiconductor's conduction band (T_n: electrons leaving the semiconductor's conduction
        // band count positive, so this pair runs from the semiconductor's side as well).
        out.J_evb = run(Band::conduction, false, EcG, &Ev, Fp, FG, t.richardson_p, -1.0,
                        out.T_p, out.d_p, 2);
        out.J_evb += run(Band::conduction, false, Ec, &EvG, Fn, FG, t.richardson_n, 1.0,
                         out.T_n, out.d_n, 2);
    }
}

DriftDiffusion::GateTunnellingCurrents DriftDiffusion::gate_tunnelling_currents(
    std::span<const double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    GateTunnellingCurrents out;
    if (gate_tunnels_.empty()) return out;
    out.ecb.assign(contact_count_, 0.0);
    out.hvb.assign(contact_count_, 0.0);
    out.evb.assign(contact_count_, 0.0);
    const std::vector<NodeDegeneracy> g = degeneracies(x);
    GateTunnelFlux flux;
    for (const GateTunnel& t : gate_tunnels_) {
        gate_tunnel_flux(t, x, g, flux);
        // T_n - T_p = c (J_ecb - J_hvb + J_evb): EVB's two directions enter T_p and T_n with the
        // signs that make both add.
        const auto c = static_cast<std::size_t>(t.contact);
        out.ecb[c] += t.coupling * flux.J_ecb;
        out.hvb[c] -= t.coupling * flux.J_hvb;
        out.evb[c] += t.coupling * flux.J_evb;
    }
    return out;
}

}  // namespace NiTCAD::assemble
