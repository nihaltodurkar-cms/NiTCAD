// The run record of a bias sweep (ARCHITECTURE.md 6.6): an identity digest of every input and the
// options as named settings.
#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

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
         {p.eps_r, p.Eg0_eV, p.varshni_alpha_eV_per_K, p.varshni_beta_K, p.Nc300, p.Nv300}) {
        d.real(v);
    }
    for (const auto& ct : {p.electron_mobility, p.hole_mobility}) {
        for (const double v : {ct.mu_min, ct.mu_max, ct.N_ref, ct.alpha, ct.T_exponent}) d.real(v);
    }
    for (const double v : {p.lifetime.tau_n0, p.lifetime.tau_p0, p.lifetime.N_ref}) d.real(v);
}

std::vector<std::pair<std::string, double>> settings(const BiasOptions& o) {
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
        {"models.doping_mobility", o.models.doping_mobility ? 1.0 : 0.0},
        {"models.srh", o.models.srh ? 1.0 : 0.0},
    };
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
        material(d, r.material.parameters());
    }
    d.integers<device::RegionId>(device.node_region());
    d.reals(device.donors());
    d.reals(device.acceptors());
    d.integer(device.contacts().size());
    for (const device::Contact& c : device.contacts()) {
        d.text(c.name);
        d.integer(static_cast<std::uint64_t>(c.kind));
        d.integers<mesh::NodeId>(c.nodes);
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
        d.reals(initial->n_cm3);
        d.reals(initial->p_cm3);
    }
    return {d.value(), named};
}

}  // namespace NiTCAD::solve
