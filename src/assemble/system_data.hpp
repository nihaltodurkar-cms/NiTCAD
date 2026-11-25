#pragma once

// The model data of a System and its state (Unit 57, DECISIONS.md GE2). Internal to the
// assemble library.
//
// SystemData is fixed at creation: the per-node, per-edge and per-contact data, the scaled tables,
// the composition and the Jacobian positions of the terms' entries. Copies of a System
// share it (shared_ptr<const SystemData>). SystemState is what changes after creation: the contact
// biases and the values they fix, the gate biases, the thermal rises, held carrier temperatures,
// the minority floor, and the tunnel paths with their layouts, the base triplets and the pattern,
// which set_paths rebuilds. Each copy owns its state, so set_bias or set_paths on one copy never
// reaches another. set_paths also moves the positions into the new pattern: it does so in a copy of
// the data (copy on write), so the data another copy shares is not changed.
//
// The terms (terms.hpp) read (const SystemData&, const SystemState&, const Evaluation&); the
// helpers they call are SystemData's members, defined in the source file of their model.

#include <complex>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/band_edges.hpp"
#include "NiTCAD/assemble/charges.hpp"
#include "NiTCAD/assemble/currents.hpp"
#include "NiTCAD/assemble/energy.hpp"
#include "NiTCAD/assemble/gate.hpp"
#include "NiTCAD/assemble/interface_edges.hpp"
#include "NiTCAD/assemble/rates.hpp"
#include "NiTCAD/assemble/sg_flux.hpp"
#include "NiTCAD/assemble/system.hpp"
#include "NiTCAD/assemble/thermal.hpp"
#include "NiTCAD/assemble/tunnelling.hpp"
#include "NiTCAD/physics/impact_ionization.hpp"
#include "NiTCAD/physics/ionization.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "NiTCAD/physics/surface_mobility.hpp"

namespace NiTCAD::assemble {

struct SystemState;

// The passkey of System's internal accessors (Unit 57 GE3): only SystemData can make one.
class SystemKey {
    SystemKey() = default;
    friend struct SystemData;
};

struct SystemData {
    using TimeStep = System::TimeStep;

    // The equation terms of the rows and their composition from the device and the models
    // (Unit 40, DECISIONS.md FE2, FE4; src/assemble/terms.hpp, composition.hpp).
    class Term;
    struct Terms;
    class Composition;
    // Fields and terms added to the drift-diffusion composition from outside its tables (Unit 57
    // GG8; composition.hpp).
    struct CompositionExtension;

    // The current density below which impact ionization fades out (system.hpp's header
    // comment). [A/cm^2]
    static constexpr double impact_current_floor = 1e-12;
    // The share of the largest rate below which a tunnel path does not matter (paths_agree), and
    // below which trace_paths drops it.
    static constexpr double matters = 1e-10;

    // A system's data and state, through System's passkey accessors.
    [[nodiscard]] static const SystemData& of(const System& s) noexcept {
        return s.data(SystemKey{});
    }
    [[nodiscard]] static const SystemState& state_of(const System& s) noexcept {
        return s.state(SystemKey{});
    }
    [[nodiscard]] static SystemState& state_of(System& s) noexcept { return s.state(SystemKey{}); }
    static void replace(System& s, std::shared_ptr<const SystemData> data) noexcept {
        s.replace_data(SystemKey{}, std::move(data));
    }

    // The model-specific queries of the topic headers (currents.hpp, thermal.hpp,
    // tunnelling.hpp, energy.hpp, charges.hpp, rates.hpp, band_edges.hpp), over the data and a
    // state; each free function of a topic header calls its own.
    [[nodiscard]] std::vector<std::pair<double, double>> edge_currents(
        const SystemState& st, std::span<const double> x) const;
    [[nodiscard]] std::vector<double> terminal_currents(const SystemState& st,
                                                        std::span<const double> x) const;
    [[nodiscard]] std::vector<double> edge_flux_terms(const SystemState& st,
                                                      std::span<const double> x) const;
    [[nodiscard]] std::vector<double> conduction_currents(const SystemState& st,
                                                          std::span<const double> x,
                                                          const TimeStep& step) const;
    [[nodiscard]] std::vector<CurrentRow> small_signal_currents(const SystemState& st,
                                                                std::span<const double> x,
                                                                std::complex<double> s) const;
    [[nodiscard]] std::vector<double> trap_occupancies(const SystemState& st,
                                                       std::span<const double> x,
                                                       const TimeStep* step) const;
    [[nodiscard]] std::vector<double> contact_charges(const SystemState& st,
                                                      std::span<const double> x,
                                                      const TimeStep* step) const;
    [[nodiscard]] std::vector<double> gate_charges(const SystemState& st,
                                                   std::span<const double> x,
                                                   const TimeStep* step) const;
    [[nodiscard]] std::vector<double> interface_trap_charges(const SystemState& st,
                                                             std::span<const double> x,
                                                             const TimeStep* step) const;
    [[nodiscard]] BandEdges band_edges(const SystemState& st, std::span<const double> x) const;
    [[nodiscard]] std::vector<MechanismRate> rates(const SystemState& st,
                                                   std::span<const double> x) const;
    [[nodiscard]] TunnelPaths trace_paths(const SystemState& st, std::span<const double> x) const;
    [[nodiscard]] bool paths_agree(const SystemState& st, const TunnelPaths& a,
                                   const TunnelPaths& b, std::span<const double> x) const;
    [[nodiscard]] std::vector<TunnelPathState> path_states(const SystemState& st,
                                                           std::span<const double> x) const;
    [[nodiscard]] GateTunnellingCurrents gate_tunnelling_currents(
        const SystemState& st, std::span<const double> x) const;
    [[nodiscard]] bool reaches_heat_sink(std::size_t node) const;
    [[nodiscard]] std::expected<void, base::Error> set_thermal_bias(
        SystemState& st, std::span<const double> temperature_K) const;
    [[nodiscard]] std::vector<double> thermal_contact_heat(const SystemState& st,
                                                           std::span<const double> x,
                                                           const TimeStep* step) const;
    [[nodiscard]] double electrical_power(const SystemState& st, std::span<const double> x) const;
    [[nodiscard]] std::vector<std::pair<double, double>> edge_joule_heat(
        const SystemState& st, std::span<const double> x) const;
    void hold_temperatures(SystemState& st, std::span<const double> x) const;

    [[nodiscard]] std::size_t node_count() const noexcept { return volume_.size(); }
    [[nodiscard]] std::size_t unknowns() const noexcept { return m_ * node_count(); }
    [[nodiscard]] std::size_t stride() const noexcept { return m_; }
    [[nodiscard]] std::size_t storage_width() const noexcept { return m_ - 1; }
    [[nodiscard]] std::size_t contact_count() const noexcept { return contact_count_; }
    [[nodiscard]] std::size_t trap_slots() const noexcept { return interfaces_.slot_count(); }
    [[nodiscard]] bool electrothermal() const noexcept { return electrothermal_; }
    [[nodiscard]] bool energy_balance() const noexcept { return energy_balance_; }
    [[nodiscard]] bool tunnelling() const noexcept {
        return tunnel_kind_ != NonlocalTunnelling::off;
    }
    [[nodiscard]] bool gate_tunnelling() const noexcept { return !gate_tunnels_.empty(); }

    struct EdgeTerm {
        std::size_t a, b;   // end nodes
        double c;           // et * scaled coupling, Poisson
        double an, ap;      // electron and hole SG edge factors (low field)
        double mu_n, mu_p;  // low-field edge mobilities (harmonic means) [cm^2/(V s)]
        double field;       // [V/cm] per unit |psi_b - psi_a|: V_T / length
        physics::CanaliParameters sat_n, sat_p;  // of node a's material
        // An edge between two materials: the ends' low-field mobilities and node b's Canali
        // parameters (field mobility takes the harmonic mean of the ends' Canali mobilities).
        bool mixed;
        double mu_n_a, mu_n_b, mu_p_a, mu_p_b;
        physics::CanaliParameters sat_n_b, sat_p_b;
        double shift_n;     // s_b - s_a + ln(n_ie,b / n_ie,a): delta_n = psi_b - psi_a + shift_n
        double shift_p;     // s_b - s_a - ln(n_ie,b / n_ie,a)
        // Thermionic emission (an interface edge with the model on): factors K and, per band,
        // ln(N_b / N_a) and N_a / N_b.
        bool thermionic;
        bool carriers;      // both ends are semiconductor nodes (else no carrier flux)
        bool charged;       // an interface edge: its Poisson flux is the interface's
        double te_kn, te_kp, te_log_nc, te_ratio_nc, te_log_nv, te_ratio_nv;
        // Positions in the Jacobian values: row a with columns of b, row b with columns of a,
        // in the order (psi, psi), (n, psi), (n, n), (p, psi), (p, p).
        std::size_t ab[5], ba[5];
        // Lombardi (Unit 46): the ends' bulk-model impurity densities [cm^-3] and the edge's
        // E_perp stencil, surface_stencil_[surface_first, surface_last); surface is false
        // without the model and on thermionic edges.
        bool surface;
        double N_a, N_b;
        std::size_t surface_first, surface_last;
    };

    // Lombardi (Unit 46): one node j of an edge's E_perp stencil, v = sum_j g_j psi_j [V/cm per
    // scaled psi], and the Jacobian positions of psi_j's column in the rows (electron, hole,
    // heat) of end a ([0]) and end b ([1]); the heat rows' only with the electrothermal model;
    // with energy balance (Unit 52) the rows electron, hole, tau_n and tau_p.
    struct SurfaceStencil {
        std::size_t node;
        double g[3];
        std::size_t pos[2][4];
    };
    // An edge's E_perp [V/cm] at a state and v, the part of its field normal to it.
    struct SurfaceField {
        double value;
        double v[3];
    };
    [[nodiscard]] SurfaceField surface_field(const EdgeTerm& e, std::span<const double> x) const;
    // d flux / d E_perp of an edge's electron and hole fluxes, with the field they were taken in.
    struct SurfaceSlope {
        SurfaceField field;
        double n, p;
    };
    // d E_perp / d psi_j of stencil entry j.
    [[nodiscard]] double surface_gradient(const SurfaceField& f, const SurfaceStencil& j) const;
    // The stencils' Jacobian positions (after the pattern is built or changed).
    void locate_surface(const linalg::SparseMatrix& pattern);
    // Lombardi's per-end low-field mobility of an edge and its slopes: the bulk mobility mu_b at
    // its impurity density N and temperature T [K] with the surface terms in E_perp.
    [[nodiscard]] static physics::SurfaceMobility surface_end(physics::Carrier carrier, double mu_b,
                                                              double N, double T, double E_perp);

    // The degeneracy factors of a node's electron and hole densities (Fermi-Dirac).
    struct NodeDegeneracy {
        physics::Degeneracy n, p;
    };

    // The storage of a semiconductor node and its derivatives d S_n / dn, d S_p / dp.
    struct NodeStorage {
        double n, d_n, p, d_p;
        double d_ln = 0.0, d_lp = 0.0;  // density gradient (Unit 62): d S / d Lambda (incomplete ionization)
    };
    [[nodiscard]] NodeStorage node_storage(std::size_t i, std::span<const double> x,
                                           std::span<const NodeDegeneracy> g) const;
    // Energy balance (Unit 52, ED9): the thermal energy the carriers of a semiconductor node store
    // per volume, W = (3/2) c theta g in units of k T_L Ns (g = F_{3/2} / F_{1/2} at c theta^(-3/2)
    // with Fermi-Dirac, 1 with Boltzmann), and its partials in the density and the rise.
    struct EnergyStorage {
        double n, n_c, n_t, p, p_c, p_t;
    };
    [[nodiscard]] EnergyStorage energy_storage(std::size_t i, std::span<const double> x) const;
    // The trap step of interface edge k, or nullptr without a step.
    [[nodiscard]] std::optional<TrapStep> trap_step(std::size_t k, const TimeStep* step) const;
    // The interface of edge k at state x.
    [[nodiscard]] InterfaceDrift interface_at(std::size_t k, std::span<const double> x,
                                              const TimeStep* step) const;

    // Impact ionization (Unit 19; see the header comment). Per node that ionizes: its carrier
    // edges (in impact_edges_) and the pseudo-inverse of sum t t^T over them (row-major, D x D);
    // per such edge: the unit vector t from its first to its second node, the field per unit
    // psi_b - psi_a and the current density per unit scaled flux, and the Jacobian positions of
    // the node's electron row with the other end's p column and of its hole row with the other
    // end's n column. Empty when the model is off.
    struct ImpactNode {
        std::size_t node;
        std::size_t first, last;
        double inverse[9];
        physics::ImpactIonizationCoefficients n, p;
        double gamma;  // the temperature factor of the node's material
        double phonon_eV;  // its optical-phonon energy (gamma at the node's T, electrothermal)
        // Energy balance (Unit 52, ED6): the node's low-field mobilities [cm^2/(V s)] and Canali
        // parameters, for E_eff(T) at the node.
        double mu_n, mu_p;
        physics::CanaliParameters sat_n, sat_p;
    };
    struct ImpactEdge {
        std::size_t edge;
        double t[3];
        double field;    // [V/cm] per unit (psi_b - psi_a): V_T / length
        double current;  // [A/cm^2] per unit scaled flux: J0 / scaled coupling area
        std::size_t n_p, p_n;
        // Density gradient (Unit 62): the node's electron row with the other end's Lambda_p
        // column and its hole row with the other end's Lambda_n column (0 with the model off).
        std::size_t n_lp = 0, p_ln = 0;
    };

    // Local band-to-band tunnelling (Unit 20): per node, its carrier edges (in field_edges_) and
    // the pseudo-inverse of sum t t^T, as ImpactNode, and its material's Kane pair.
    struct TunnelNode {
        std::size_t node;
        std::size_t first, last;
        double inverse[9];
        double A, B;
    };
    struct FieldEdge {
        std::size_t edge;
        double t[3];
        double field;  // [V/cm] per unit (psi_b - psi_a)
    };

    // Nonlocal tunnelling (Unit 20): per frozen path, the nodes whose psi its rows read (columns),
    // the nodes its electrons go to (deposit), the start's hole-row positions per column and the
    // deposit rows' positions (deposit-major), and the forward length before each segment [cm].
    struct PathLayout {
        std::vector<std::size_t> columns, deposit;
        std::vector<std::size_t> hole, electron;
        std::vector<double> before_cm;
    };
    // One path's live terms at x: its pairs P (scaled) and dP / dpsi per column, the share of each
    // deposit node and its partials (deposit-major), and the diagnostics.
    struct PathTerms {
        double pairs = 0.0;
        std::vector<double> d_pairs;
        std::vector<double> share, d_share;
        TunnelPathState state{};
    };
    // The columns, deposit nodes and lengths of a path (no Jacobian positions).
    [[nodiscard]] PathLayout layout(const TunnelPath& p) const;
    void path_terms(const TunnelPath& p, const PathLayout& l, std::span<const double> x,
                    PathTerms& t) const;
    // The Jacobian positions every term uses, after the pattern changed (set_paths).
    void locate_positions(const linalg::SparseMatrix& pattern);

    // Gate tunnelling (Unit 51; gate_tunnels.cpp, DECISIONS.md GD1-GD8). One tunnel per
    // semiconductor node of a tunnelling gate (lumped) or per semiconductor-to-insulator edge whose
    // straight grid line reaches a tunnelling electrode. The barrier's points are vacuum levels
    // depth - V_T psi at psi = w psi_a + (1 - w) psi_b of two local columns (a == b: one node), or
    // (lumped == true, the last point) the gate's own, -V_G + phi_m. Local columns: 0, 1, 2 the
    // node's psi, n, p; 3 + k the psi of nodes[k]; the bias of a lumped gate after them.
    struct TunnelPoint {
        std::size_t a, b;
        double w;
    };
    struct TunnelPiece {
        std::size_t p0, p1;   // points
        double length_cm;
        double affinity_eV, gap_eV, mass_n, mass_p, mass_evb;
    };
    struct GateTunnel {
        std::size_t node;
        std::int32_t contact;
        std::uint8_t processes;
        bool lumped;
        device::GateElectrode electrode;
        double work_function_eV;  // phi_m of the electrode
        double gate_gap_eV;       // a polysilicon electrode's E_g(T)
        double coupling;          // scaled particle flow per A/cm^2: A / (q Ns L_D^(D-2) D0)
        double richardson_n, richardson_p;  // A* of the node's bands [A/(cm^2 K^2)]
        double depth_eV;          // the reference material's intrinsic-level depth
        std::vector<std::size_t> nodes;  // path nodes, the electrode node last (empty: lumped)
        std::vector<TunnelPoint> points;
        std::vector<TunnelPiece> pieces;
        // Jacobian positions of the node's electron and hole rows with the psi of each path node.
        std::vector<std::size_t> pos_n, pos_p;
    };
    // A tunnel's row terms at x: the scaled sinks T_n (electron row: -T_n) and T_p (hole row:
    // +T_p), the gate's terminal current T_n - T_p, and their partials per local column.
    struct GateTunnelFlux {
        double T_n = 0.0, T_p = 0.0;
        std::vector<double> d_n, d_p;
        double J_ecb = 0.0, J_hvb = 0.0, J_evb = 0.0;  // [A/cm^2], diagnostics
    };
    [[nodiscard]] std::optional<base::Error> make_gate_tunnels(const device::Device& device,
                                                               const Scaling& scaling);
    void locate_gate_tunnels(const linalg::SparseMatrix& pattern);
    void gate_tunnel_flux(const SystemState& st, const GateTunnel& t, std::span<const double> x,
                          std::span<const NodeDegeneracy> g, GateTunnelFlux& out) const;
    // Every tunnel's flux at x, out[k] gate_tunnels_[k]'s, the tunnels shared among threads (the
    // first threading inside assembly, GR11). Each flux is computed whole by one thread, so the
    // values do not depend on the thread count; callers add them in tunnel order.
    void gate_tunnel_fluxes(const SystemState& st, std::span<const double> x,
                            std::span<const NodeDegeneracy> g,
                            std::vector<GateTunnelFlux>& out) const;
    // The same, started on worker threads when an assembly starts and collected by the gate
    // tunnelling term, so the other terms add theirs meanwhile (terms.hpp).
    class GateTunnelJob;

    // Per node at state x under Fermi-Dirac statistics; empty under Boltzmann.
    [[nodiscard]] std::vector<NodeDegeneracy> degeneracies(std::span<const double> x) const;

    // The electron and hole fluxes of an edge at state x, with field-dependent edge factors when
    // the model is on and the degeneracy terms g (degeneracies(x)) when it is not empty. With
    // Lombardi on the edge (e.surface) the fluxes' partials exclude E_perp's; slope, when not
    // nullptr, receives d flux / d E_perp (0 without the model).
    [[nodiscard]] std::pair<EdgeFlux, EdgeFlux> edge_fluxes(
        const EdgeTerm& e, std::span<const double> x, std::span<const NodeDegeneracy> g,
        SurfaceSlope* slope = nullptr) const;

    // Electrothermal (Unit 23; electrothermal.cpp). Per node, constant: the material (index into
    // thermal_materials_, or -1 on an insulator), its thermal data, the impurity densities the
    // mobility and band-gap narrowing read, the band constants c_n, c_p and Eg(T0), the
    // thermopower exponents, ln Nc0 and ln Nv0 (scaled), the thermal contact (-1: none), the
    // R_th coefficient h (0: none) and the heat capacity C = rho c V L_D^2 / (t0 kappa_ref) of the
    // heat row (times the rate and the change of tau).
    struct ThermalNode {
        std::int32_t material;
        physics::ThermalParameters thermal;
        double impurity, impurity_mobility, narrowing_eV;
        double c_n, c_p, gap0_eV, r_n, r_p, log_nc, log_nv;
        std::int32_t contact;
        double h;
        double capacity;
        // A p-polysilicon work function, the only one that moves with T (chi + Eg(T)): 1 on a
        // node of a gate with one (its electrode potential gains the node's own g(theta)), 2 on a
        // node of an electrode of one (it gains polysilicon's, silicon's -(Eg(T) - Eg(T0)) / k T0).
        std::int32_t ppoly;
        // Unit 66: N_D, N_A for the Klaassen mobility (0 with doping_mobility off) [cm^-3].
        double donors_mobility = 0.0, acceptors_mobility = 0.0;
    };
    // Per node at a state: theta and tau = theta - 1 (the unknown), n_ie / Ns, the mobilities
    // [cm^2/(V s)], g (the valence band's shift), kappa [W/(cm K)] and the Kirchhoff transform u,
    // each with d / d theta.
    // mass_n, mass_p: ln M(theta) = ln((m_d / m0)^(3/2) at T / at T0) of the conduction and
    // valence bands (Unit 35; 0 without a density-of-states mass fit), d_mass_* d / d theta.
    struct ThermalState {
        double theta, tau, nie, d_nie, mu_n, d_mu_n, mu_p, d_mu_p, gap, d_gap, kappa, d_kappa, u;
        double mass_n, d_mass_n, mass_p, d_mass_p;
    };
    // Per edge, constant: V_T / D0 times the coupling (g_e), L_c, L_v, the mean thermopower
    // exponents, whether conduction takes the Kirchhoff form (both ends in one material), g_e /
    // kappa_ref, and the Jacobian positions of the full 4 x 4 blocks (row a with columns of b,
    // row b with columns of a; row-major).
    struct ThermalEdge {
        double coupling, log_nc, log_nv, r_n, r_p;
        bool kirchhoff;
        double conduction;
        std::size_t ab[16], ba[16];
    };
    // A flux and its partials with respect to psi, the density and theta of each end (a, b).
    struct ThermalFlux {
        double flux;
        double d_psi[2], d_c[2], d_theta[2];
        double delta, theta_e, w[2];  // for the Joule part
    };
    [[nodiscard]] std::vector<ThermalState> thermal_states(std::span<const double> x) const;
    // A carrier's band reduced energy x = ln(c / N(theta)) - ln gamma (Fermi-Dirac; ln gamma 0
    // under Boltzmann statistics) at its node's T, with ln gamma and the Fermi-Dirac thermal
    // diffusion factor's excess over the Boltzmann 1 + r (r = -1/2): h - 1/2, each with its
    // partials in the density c and theta. log_N: ln N0 (scaled) of the band.
    struct Level {
        double x, dx_dc, dx_dt;
        double log_gamma, dg_dc, dg_dt;
        double extra, dh_dc, dh_dt;
    };
    // mass, d_mass: the band's ln M(theta) and its theta derivative (ThermalState).
    [[nodiscard]] Level level(double density, double theta, double nie, double log_N, double mass,
                              double d_mass) const;
    struct NodeLevels {
        Level n, p;
    };
    // Per node at x under Fermi-Dirac statistics; empty under Boltzmann.
    [[nodiscard]] std::vector<NodeLevels> thermal_levels(std::span<const double> x,
                                                         std::span<const ThermalState> t) const;
    // An edge's electron and hole fluxes with the electrothermal model; slope as edge_fluxes'.
    [[nodiscard]] std::pair<ThermalFlux, ThermalFlux> thermal_fluxes(
        std::size_t edge, std::span<const double> x, std::span<const ThermalState> t,
        std::span<const NodeLevels> g, SurfaceSlope* slope = nullptr) const;
    // The heat row's conduction term along an edge, from a to b, and its theta partials.
    struct Conduction {
        double flux, d_a, d_b;
    };
    [[nodiscard]] Conduction conduction(std::size_t edge, std::span<const ThermalState> t) const;
    // An ohmic contact node's Dirichlet values at its T and their theta derivatives.
    struct ThermalOhmic {
        double psi, n, p, d_psi, d_n, d_p;
    };
    [[nodiscard]] ThermalOhmic thermal_ohmic(const SystemState& st, std::size_t node,
                                             const ThermalState& t) const;
    // The shift of a gate or electrode node's electrode potential from its value at T0 (ppoly),
    // and its theta derivative; 0 for every other node.
    struct Shift {
        double value, d_theta;
    };
    [[nodiscard]] Shift electrode_shift(std::size_t node, const ThermalState& t) const;
    // Adds G shift of every p-polysilicon gate node to its gate's charge (gates_.charges is at T0).
    void add_gate_shift(const SystemState& st, std::span<const double> x,
                        std::span<double> charge) const;
    // Creation of the thermal data; errors as create's (see electrothermal.cpp).
    // edge_geometry: each mesh edge's scaled coupling (coupling area / length / L_D^(D-2)).
    [[nodiscard]] std::optional<base::Error> make_thermal(SystemState& st,
                                                          const device::Device& device,
                                                          const Scaling& scaling,
                                                          const PhysicsModels& models,
                                                          std::span<const double> edge_geometry);

    // The unknowns (psi, n, p, with the electrothermal model theta last) and the rows' terms in
    // assembly order (composition.hpp; DECISIONS.md FE3, FE4). Immutable, shared by copies.
    std::shared_ptr<const Composition> composition_;
    std::size_t m_ = 3;  // unknowns per node: composition_->dofs().block()
    // The built-in fields of each node block (m_ less a CompositionExtension's fields, which follow
    // them), and whether the heat equation is in (models.electrothermal).
    std::size_t built_in_ = 3;
    bool electrothermal_ = false;
    std::vector<physics::Semiconductor> thermal_materials_;
    std::vector<ThermalNode> thermal_nodes_;
    std::vector<ThermalEdge> thermal_edges_;
    // A thermal contact's kind (true: isothermal) and its nodes' indices into its list.
    std::vector<char> thermal_isothermal_;
    double T0_ = 0.0, heat_scale_ = 0.0, kappa_ref_ = 0.0;
    bool thermal_bgn_ = false;
    // Units 34, 35: per thermal material, whether its hole mobility is Masetti's.
    std::vector<char> masetti_holes_;
    bool masetti_electrons_ = false;  // Unit 38: the electron mobility is Masetti's
    // Unit 66: the electron / hole mobility is Klaassen's at neutral equilibrium densities.
    bool klaassen_electrons_ = false, klaassen_holes_ = false;
    bool sinks_complete_ = false;
    std::vector<char> reaches_sink_;  // per node (reaches_heat_sink)
    // One node's thermal state at tau (NaN entries unless tau is finite and theta positive).
    [[nodiscard]] ThermalState thermal_state(std::size_t node, double tau) const;
    // s times the Jacobian of the electrothermal storage terms (assemble_heat with a step) about a
    // steady state, added to the small-signal matrix's values.
    void add_thermal_storage(std::span<const double> x, std::complex<double> s,
                             std::span<std::complex<double>> values) const;
    // The electrothermal rows (assemble with the model on); with sinks false, without the thermal
    // contacts' terms, so every heat row is the node's heat balance alone.
    void assemble_heat(const SystemState& st, std::span<const double> x, std::span<double> f,
                       std::span<double> values, const TimeStep* step, bool sinks) const;

    std::vector<double> volume_, doping_, n_ie_, tau_n_, tau_p_, auger_n_, auger_p_;
    std::vector<double> log_dos_n_, log_dos_p_;  // ln(Nc / n_ie), ln(Nv / n_ie)
    std::vector<double> band_shift_;             // s
    std::vector<double> donors_, acceptors_;     // N_D / Ns, N_A / Ns
    std::vector<physics::DopantLevels> levels_;  // dopant levels in units of kT
    std::vector<double> radiative_;              // B [cm^3/s]
    std::vector<std::int32_t> contact_;   // ohmic contact index per node, or -1
    // Unit 48: the Schottky contact index per node, or -1; eta_B, K_n and K_p per node
    // (scaled_device.hpp); and the contact whose terminal current and charge a node's edges count
    // toward (ohmic or Schottky), or -1.
    std::vector<std::int32_t> schottky_;
    std::vector<double> schottky_eta_, schottky_kn_, schottky_kp_;
    std::vector<std::int32_t> owner_;
    std::vector<char> insulator_;         // 1 on insulator nodes
    std::vector<std::int32_t> electrode_; // electrode contact index per node, or -1
    std::vector<double> electrode_potential_;  // psi_E at zero bias
    InterfaceEdges interfaces_;
    // Per interface edge, the Jacobian positions of the insulator node's Poisson row with the
    // semiconductor node's n and p columns.
    std::vector<std::pair<std::size_t, std::size_t>> interface_np_;

    [[nodiscard]] InterfaceStatistics statistics(std::size_t node) const noexcept;
    std::vector<device::ContactKind> kinds_;  // per contact
    std::vector<EdgeTerm> edges_;
    std::vector<std::size_t> block_;      // per node, 9 positions of its 3x3 diagonal block
    std::vector<ImpactNode> impact_nodes_;
    std::vector<ImpactEdge> impact_edges_;
    std::vector<TunnelNode> tunnel_nodes_;
    std::vector<FieldEdge> field_edges_;
    // Nonlocal tunnelling: the model, the device's cells, per node the longest path that may start
    // there (0: none) [cm], the material's gap at T [eV] and tunnelling masses [kg] (direct_wkb),
    // and the Kane pair (kane). The frozen paths are SystemState's.
    NonlocalTunnelling tunnel_kind_ = NonlocalTunnelling::off;
    std::optional<mesh::TensorCells> cells_;
    std::vector<double> tunnel_length_, tunnel_gap_, tunnel_mc_, tunnel_mv_, tunnel_A_, tunnel_B_;
    std::vector<GateTunnel> gate_tunnels_;  // Unit 51; empty without a tunnelling contact
    double temperature_K_ = 0.0;
    bool generation_ = false;  // any generation model on
    int dimension_ = 1;
    double generation_scale_ = 0.0;        // 1 / R0: [cm^-3 s^-1] to the rows' scaled rate
    const ManufacturedSources* manufactured_ = nullptr;  // models.manufactured (TVQ3), not owned
    double impact_resolution_ = 0.0;       // models.impact_current_resolution
    std::size_t contact_count_ = 0;
    double rate_scale_ = 0.0;             // Ns / R0: scaled SRH call to R / R0
    double Ns_ = 0.0;
    double V_T_ = 0.0;
    bool srh_ = true;
    bool auger_ = true;
    bool field_mobility_ = false;
    // Lombardi (Unit 46): on, the device's temperature [K] (isothermal) and the edges' stencils.
    bool surface_ = false;
    double surface_T_ = 0.0;
    std::vector<SurfaceStencil> surface_stencil_;
    bool fermi_dirac_ = false;
    bool ionization_ = false;
    bool radiative_on_ = true;

    // Energy balance (Unit 52; energy_balance.cpp). Per node the energy relaxation times [s] (0
    // off the semiconductor); per mesh edge (parallel to edges_) those of its ends and the
    // positions of its full 5 x 5 blocks, row a with the columns of b (ab) and row b with a's (ba).
    bool energy_balance_ = false;
    std::vector<double> relax_n_, relax_p_;
    struct EnergyEdge {
        std::size_t ab[25], ba[25];
    };
    std::vector<EnergyEdge> energy_edges_;
    // One carrier's particle flux on an edge (the continuity rows' quantity, edge_fluxes' sign)
    // and its energy flux into end a [k T_L units], each with partials in psi, the density and
    // the temperature rise of ends a ([0]) and b ([1]); delta the edge's driving term (the band
    // edge's drop from b to a in k T_L, sign as the flux's), for the Joule term -delta flux / 2.
    // Lombardi (surface edges): d_perp and e_perp, the fluxes' partials in the edge's E_perp.
    struct EnergyFlux {
        double flux, d_psi[2], d_c[2], d_t[2];
        double energy, e_psi[2], e_c[2], e_t[2];
        double delta;
        double d_perp, e_perp;
    };
    [[nodiscard]] std::pair<EnergyFlux, EnergyFlux> energy_fluxes(std::size_t k,
                                                                  std::span<const double> x) const;
    // Creation of the energy-balance data after the pattern; errors as create's.
    [[nodiscard]] std::optional<base::Error> make_energy_balance(SystemState& st,
                                                                 const device::Device& device,
                                                                 const PhysicsModels& models);

    // Density gradient (Unit 62; DECISIONS.md QP2-QP5, QP8, QP13; density_gradient.cpp). The
    // quantum potentials Lambda_n / kT and Lambda_p / kT are unknowns after psi, n, p (stride 5) or,
    // at thermal equilibrium, after psi alone (stride 3); every branch reads this flag and the
    // slots below, never the stride. beta / L_D^2 per carrier, the barrier terms per node
    // (scaled_device.hpp), and per mesh edge (parallel to edges_) its geometry g = coupling area /
    // length / L_D^(D-2) and, on carrier edges, the Jacobian positions of row a with the columns of
    // b (ab) and of row b with a's (ba): drift-diffusion (Lambda_n, n), (Lambda_p, p),
    // (n, Lambda_n), (p, Lambda_p); equilibrium (Lambda_n, psi), (Lambda_n, Lambda_n),
    // (Lambda_p, psi), (Lambda_p, Lambda_p). Per interface edge, the insulator node's Poisson row
    // with the semiconductor node's Lambda_n and Lambda_p columns (drift-diffusion).
    bool density_gradient_ = false;
    bool dg_impact_ = false;  // impact ionization on (its cross-carrier entries)
    double dg_beta_n_ = 0.0, dg_beta_p_ = 0.0;
    std::vector<double> dg_barrier_n_, dg_barrier_p_, dg_geometry_;
    struct DgEdge {
        std::size_t ab[4], ba[4];
    };
    std::vector<DgEdge> dg_edges_;
    std::vector<std::pair<std::size_t, std::size_t>> interface_lambda_;
    [[nodiscard]] bool density_gradient() const noexcept { return density_gradient_; }
    // The slots of Lambda_n and Lambda_p in a node block.
    [[nodiscard]] std::size_t lambda_n_slot() const noexcept { return equilibrium_ ? 1 : 3; }
    [[nodiscard]] std::size_t lambda_p_slot() const noexcept { return equilibrium_ ? 2 : 4; }
    // The pattern entries beyond every node's full block and the edges' stencil (both
    // formulations; with dg_impact_ the cross-carrier entries of impact ionization), and their
    // positions (after the impact edges are made).
    void density_gradient_pattern(std::vector<linalg::Triplet>& t) const;
    void locate_density_gradient(const linalg::SparseMatrix& pattern);
    // Node i's statistics with its quantum potentials at x (InterfaceStatistics; drift-diffusion).
    [[nodiscard]] InterfaceStatistics statistics(std::size_t node,
                                                 std::span<const double> x) const noexcept;

    // Thermal equilibrium (Unit 57 GE5; equilibrium.cpp, equilibrium.hpp): psi alone (m_ = 1),
    // the terms of equilibrium_table. Of the data only the scaled device's tables are filled;
    // the edges carry c, charged, carriers and the positions ab[0], ba[0]; the state the gates
    // and psi0_ (the Dirichlet values: an ohmic contact's at zero bias, a Schottky contact's
    // psi_S at zero bias, an electrode's at its bias).
    bool equilibrium_ = false;
    struct Made;
    [[nodiscard]] static std::expected<Made, base::Error> make_equilibrium(
        const device::Device& device, const Scaling& scaling, const PhysicsModels& models);
    // The coupled drift-diffusion system (System::create's, drift_diffusion.cpp), before its
    // zero bias is set; with an extension, its fields follow the built-in ones in every node
    // block and its terms follow the table's (composition.hpp).
    [[nodiscard]] static std::expected<Made, base::Error> make_drift_diffusion(
        const device::Device& device, const Scaling& scaling, const PhysicsModels& models,
        const CompositionExtension* extension);
    // System::create with an extension (Unit 57 GG8), at zero bias. Precondition
    // (NITCAD_EXPECTS): the isothermal drift-diffusion models (neither electrothermal nor energy
    // balance), at least one extension field, each a node unknown.
    [[nodiscard]] static std::expected<System, base::Error> create_extended(
        const device::Device& device, const Scaling& scaling, const PhysicsModels& models,
        const CompositionExtension& extension);
    [[nodiscard]] std::expected<void, base::Error> set_equilibrium_bias(
        SystemState& st, std::span<const double> bias_V) const;
    // The carriers slaved to psi (0 on insulator nodes).
    void equilibrium_carriers(std::span<const double> psi, std::span<double> n,
                              std::span<double> p) const;
    [[nodiscard]] std::vector<double> charge_neutral_potential(const SystemState& st) const;
    [[nodiscard]] BandEdges equilibrium_band_edges(std::span<const double> psi) const;
    [[nodiscard]] std::vector<double> equilibrium_gate_charges(const SystemState& st,
                                                               std::span<const double> psi) const;
    [[nodiscard]] std::vector<double> equilibrium_trap_charges(std::span<const double> psi) const;
};

struct SystemState {
    GateNodes gates_;
    // Dirichlet values per node (contact nodes only); on a Schottky node psi_S, n_B and p_B.
    std::vector<double> psi0_, n0_, p0_;
    std::vector<double> metal_;          // per node: V / V_T on an ohmic contact node, else 0
    std::vector<double> contact_bias_;   // per contact: V / V_T
    std::vector<double> thermal_rise_;   // per thermal contact: (T_c - T0) / T0
    std::vector<double> held_;  // hold_temperatures: tau_n, tau_p per node (empty: not held)
    bool minority_floor_ = false;  // Unit 50 RD2, set_minority_floor
    // Nonlocal tunnelling: the frozen paths, their layouts, and the pattern's triplets without the
    // terms' own entries (Term::pattern).
    TunnelPaths paths_;
    std::vector<SystemData::PathLayout> layouts_;
    std::vector<linalg::Triplet> base_triplets_;
    linalg::SparseMatrix pattern_;
};

// A system's data and state as creation makes them (make_equilibrium).
struct SystemData::Made {
    std::shared_ptr<SystemData> data;
    SystemState state;
};

}  // namespace NiTCAD::assemble
