// The run record of a bias sweep (ARCHITECTURE.md 6.6): an identity digest of every input and the
// options as named settings.
#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/interface_traps.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"

namespace NiTCAD::solve {

namespace {

// 64-bit FNV-1a over a sequence of tagged values. Doubles are hashed by their bit pattern, so
// equal inputs give equal digests and -0.0 differs from 0.0.
class Digest {
public:
    void bytes(const void* data, std::size_t size) {
        const auto* p = static_cast<const unsigned char*>(data);
        for (std::size_t i = 0; i < size; ++i) {
            h_ ^= p[i];
            h_ *= 0x100000001b3ull;
        }
    }
    void integer(std::uint64_t v) { bytes(&v, sizeof v); }
    void real(double v) { integer(std::bit_cast<std::uint64_t>(v)); }
    void text(std::string_view s) {
        integer(s.size());
        bytes(s.data(), s.size());
    }
    void reals(std::span<const double> v) {
        integer(v.size());
        for (const double x : v) real(x);
    }
    template <class Int>
    void integers(std::span<const Int> v) {
        integer(v.size());
        for (const Int x : v) integer(static_cast<std::uint64_t>(static_cast<std::int64_t>(x)));
    }
    [[nodiscard]] std::uint64_t value() const noexcept { return h_; }

private:
    std::uint64_t h_ = 0xcbf29ce484222325ull;
};

void material(Digest& d, const physics::SemiconductorParameters& p) {
    for (const double v :
         {p.eps_r, p.Eg0_eV, p.varshni_alpha_eV_per_K, p.varshni_beta_K, p.Nc300, p.Nv300,
          p.electron_affinity_eV}) {
        d.real(v);
    }
    for (const auto& ct : {p.electron_mobility, p.hole_mobility}) {
        for (const double v : {ct.mu_min, ct.mu_max, ct.N_ref, ct.alpha, ct.T_exponent}) d.real(v);
    }
    for (const double v : {p.lifetime.tau_n0, p.lifetime.tau_p0, p.lifetime.N_ref}) d.real(v);
    for (const double v : {p.auger.Cn, p.auger.Cp}) d.real(v);
    for (const double v : {p.bandgap_narrowing.E0_eV, p.bandgap_narrowing.N0}) d.real(v);
    for (const auto& c : {p.electron_saturation, p.hole_saturation}) {
        d.real(c.v_sat_cm_s);
        d.real(c.beta);
    }
    d.real(p.radiative_cm3_s);
    for (const double v : {p.ionization.donor_eV, p.ionization.acceptor_eV,
                           p.ionization.donor_degeneracy, p.ionization.acceptor_degeneracy,
                           p.richardson.electron, p.richardson.hole}) {
        d.real(v);
    }
}

// The options the sweep reads: the quasi-static sweep has no continuity equations, so the
// mobility, SRH, Auger, field-mobility and radiative switches do not apply to it and are left out
// (of the digest too).
std::vector<std::pair<std::string, double>> settings(const BiasOptions& o) {
    const bool transport = o.equations == Equations::drift_diffusion;
    std::vector<std::pair<std::string, double>> s{
        {"newton.max_iterations", static_cast<double>(o.newton.max_iterations)},
        {"newton.tol_update", o.newton.tol_update},
        {"newton.max_update", o.newton.max_update},
        {"linear.backend", static_cast<double>(o.linear.backend)},
        {"linear.threads", static_cast<double>(o.linear.threads)},
        {"linear.equilibrate", o.linear.equilibrate ? 1.0 : 0.0},
        {"linear.max_backward_error", o.linear.max_backward_error},
        {"linear.max_refinement_steps", static_cast<double>(o.linear.max_refinement_steps)},
        {"linear.min_pivot_ratio", o.linear.min_pivot_ratio},
        {"models.bgn", o.models.bgn ? 1.0 : 0.0},
        {"models.fermi_dirac", o.models.fermi_dirac ? 1.0 : 0.0},
        {"models.incomplete_ionization", o.models.incomplete_ionization ? 1.0 : 0.0},
        {"equations", static_cast<double>(o.equations)},
    };
    if (transport) {
        s.emplace_back("models.doping_mobility", o.models.doping_mobility ? 1.0 : 0.0);
        s.emplace_back("models.srh", o.models.srh ? 1.0 : 0.0);
        s.emplace_back("models.auger", o.models.auger ? 1.0 : 0.0);
        s.emplace_back("models.field_mobility", o.models.field_mobility ? 1.0 : 0.0);
        s.emplace_back("models.radiative", o.models.radiative ? 1.0 : 0.0);
    }
    if (o.Ns_override) s.emplace_back("scaling.Ns_override", *o.Ns_override);
    return s;
}

}  // namespace

results::RunRecord make_run_record(const device::Device& device, const BiasOptions& options,
                                   std::span<const std::vector<double>> points,
                                   const results::NodeFields* initial) {
    Digest d;
    const mesh::Mesh& m = device.mesh();
    d.integer(static_cast<std::uint64_t>(m.dimension()));
    d.integer(m.node_count());
    for (const mesh::Point& p : m.points()) d.reals(p);
    d.reals(m.volumes());
    d.integer(m.edges().size());
    for (const mesh::Edge& e : m.edges()) {
        d.integer(static_cast<std::uint64_t>(e.first));
        d.integer(static_cast<std::uint64_t>(e.second));
        d.real(e.length);
        d.real(e.coupling_area);
    }
    d.integer(m.boundary().size());
    for (const mesh::BoundaryPatch& b : m.boundary()) {
        d.text(b.name);
        d.integers<mesh::NodeId>(b.nodes);
        d.reals(b.areas);
    }
    d.real(device.temperature_K());
    d.integer(device.regions().size());
    for (const device::Region& r : device.regions()) {
        d.text(r.name);
        if (device::is_insulator(r)) {  // tagged; a semiconductor region hashes as before
            d.text("insulator");
            d.real(std::get<physics::Insulator>(r.material).parameters().eps_r);
            continue;
        }
        material(d, std::get<physics::Semiconductor>(r.material).parameters());
    }
    d.integers<device::RegionId>(device.node_region());
    d.reals(device.donors());
    d.reals(device.acceptors());
    d.integer(device.contacts().size());
    for (const device::Contact& c : device.contacts()) {
        d.text(c.name);
        d.integer(static_cast<std::uint64_t>(c.kind));
        d.integers<mesh::NodeId>(c.nodes);
        if (c.kind == device::ContactKind::gate) {
            d.text(c.gate.boundary);
            d.integer(static_cast<std::uint64_t>(c.gate.electrode));
            for (const double v : {c.gate.oxide_thickness_cm, c.gate.oxide_relative_permittivity,
                                   c.gate.fixed_charge_cm2}) {
                d.real(v);
            }
            // A polysilicon electrode takes its work function from the semiconductor.
            if (c.gate.electrode == device::GateElectrode::metal) d.real(c.gate.work_function_eV);
        }
        if (c.kind == device::ContactKind::electrode) {
            d.integer(static_cast<std::uint64_t>(c.electrode.kind));
            if (c.electrode.kind == device::GateElectrode::metal) {
                d.real(c.electrode.work_function_eV);
            }
        }
    }
    // Interface charge and traps act on Poisson (both equation sets), recombination on the
    // continuity equations; hashed only when an interface has them, so a device without keeps its
    // digest.
    for (std::size_t k = 0; k < device.interfaces().size(); ++k) {
        const device::Interface& f = device.interfaces()[k];
        if (!device::has_interface_charge_or_recombination(f)) continue;
        d.integer(k);
        d.real(f.fixed_charge_cm2);
        const physics::InterfaceTraps& t = f.traps;
        d.integer(t.levels.size());
        for (const physics::TrapLevel& l : t.levels) {
            d.integer(static_cast<std::uint64_t>(l.type));
            for (const double v : {l.density_cm2, l.energy_eV, l.sigma_n_cm2, l.sigma_p_cm2}) {
                d.real(v);
            }
        }
        d.integer(t.bands.size());
        for (const physics::TrapBand& b : t.bands) {
            d.integer(static_cast<std::uint64_t>(b.type));
            for (const double v : {b.density_cm2_eV, b.energy_low_eV, b.energy_high_eV,
                                   b.sigma_n_cm2, b.sigma_p_cm2}) {
                d.real(v);
            }
        }
        d.real(t.thermal_velocity_n_cm_s);
        d.real(t.thermal_velocity_p_cm_s);
        if (options.equations == Equations::drift_diffusion) {
            d.real(f.recombination_velocity_n_cm_s);
            d.real(f.recombination_velocity_p_cm_s);
        }
    }
    // Interface transport acts on the continuity equations only.
    if (options.equations == Equations::drift_diffusion) {
        d.integer(device.interfaces().size());
        for (const device::Interface& f : device.interfaces()) {
            d.text(f.region_a);
            d.text(f.region_b);
            d.integer(static_cast<std::uint64_t>(f.transport));
        }
    }
    const auto named = settings(options);
    d.integer(named.size());
    for (const auto& [name, value] : named) {
        d.text(name);
        d.real(value);
    }
    d.integer(points.size());
    for (const auto& p : points) d.reals(p);
    d.integer(initial != nullptr ? 1 : 0);
    if (initial != nullptr) {
        d.reals(initial->potential_V);
        if (options.equations == Equations::drift_diffusion) {  // the quasi-static sweep reads
            d.reals(initial->n_cm3);                            // only the potential
            d.reals(initial->p_cm3);
        }
    }
    return {d.value(), named};
}

}  // namespace NiTCAD::solve
