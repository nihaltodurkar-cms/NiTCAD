#include "NiTCAD/deck/build.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <variant>

#include "NiTCAD/device/doping_profiles.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/mesh/complex.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "values.hpp"
#include "waveforms.hpp"

// Unwraps an expected into `var`, or returns its error from the enclosing function.
#define NITCAD_DECK_TRY(var, expr)                                                \
    auto var##_expected = (expr);                                                 \
    if (!var##_expected) return std::unexpected(std::move(var##_expected.error())); \
    auto var = std::move(*var##_expected)

// Returns the error of an expected<void> from the enclosing function.
#define NITCAD_DECK_CHECK(expr) \
    if (auto check_ = (expr); !check_) return std::unexpected(std::move(check_.error()))

namespace NiTCAD::deck {

namespace {

using namespace detail;

// --- material parameters by name (DK7) ---------------------------------------------------------

#define NITCAD_THERMAL_PARAMETERS(X)                                                           \
    X(thermal.conductivity_W_cmK) X(thermal.conductivity_exponent) X(thermal.heat_capacity_J_cm3K) \
    X(thermal.thermopower_exponent_n) X(thermal.thermopower_exponent_p)

#define NITCAD_SEMICONDUCTOR_PARAMETERS(X)                                                     \
    X(eps_r) X(Eg0_eV) X(varshni_alpha_eV_per_K) X(varshni_beta_K) X(Nc300) X(Nv300)           \
    X(electron_affinity_eV) X(electron_mobility.mu_min) X(electron_mobility.mu_max)            \
    X(electron_mobility.N_ref) X(electron_mobility.alpha) X(electron_mobility.T_exponent)      \
    X(hole_mobility.mu_min) X(hole_mobility.mu_max) X(hole_mobility.N_ref)                     \
    X(hole_mobility.alpha) X(hole_mobility.T_exponent) X(lifetime.tau_n0) X(lifetime.tau_p0)   \
    X(lifetime.N_ref) X(auger.Cn) X(auger.Cp) X(bandgap_narrowing.E0_eV)                       \
    X(bandgap_narrowing.N0) X(electron_saturation.v_sat_cm_s) X(electron_saturation.beta)      \
    X(hole_saturation.v_sat_cm_s) X(hole_saturation.beta) X(radiative_cm3_s)                   \
    X(ionization.donor_eV) X(ionization.acceptor_eV) X(ionization.donor_degeneracy)            \
    X(ionization.acceptor_degeneracy) X(richardson.electron) X(richardson.hole)                \
    X(impact_ionization.electron.A_low_per_cm) X(impact_ionization.electron.B_low_V_per_cm)    \
    X(impact_ionization.electron.A_high_per_cm) X(impact_ionization.electron.B_high_V_per_cm)  \
    X(impact_ionization.electron.switch_V_per_cm) X(impact_ionization.hole.A_low_per_cm)       \
    X(impact_ionization.hole.B_low_V_per_cm) X(impact_ionization.hole.A_high_per_cm)           \
    X(impact_ionization.hole.B_high_V_per_cm) X(impact_ionization.hole.switch_V_per_cm)        \
    X(impact_ionization.phonon_energy_eV) X(band_to_band.A_per_cm3_s)                          \
    X(band_to_band.B_V_per_cm) X(band_to_band.electron_mass) X(band_to_band.hole_mass)         \
    X(passler_gap.alpha_eV_per_K) X(passler_gap.Theta_K) X(passler_gap.Delta)                  \
    X(energy_relaxation.tau_n_s) X(energy_relaxation.tau_p_s)                                  \
    NITCAD_THERMAL_PARAMETERS(X)

double* semiconductor_parameter(physics::SemiconductorParameters& p, std::string_view name) {
#define NITCAD_X(path) \
    if (name == #path) return &p.path;
    NITCAD_SEMICONDUCTOR_PARAMETERS(NITCAD_X)
#undef NITCAD_X
    return nullptr;
}

double* insulator_parameter(physics::InsulatorParameters& p, std::string_view name) {
#define NITCAD_X(path) \
    if (name == #path) return &p.path;
    NITCAD_X(eps_r)
    NITCAD_THERMAL_PARAMETERS(NITCAD_X)
    // Unit 51 (GD5, GD10): gate-tunnelling band data.
    NITCAD_X(electron_affinity_eV)
    NITCAD_X(band_gap_eV)
    NITCAD_X(tunnel_mass_n)
    NITCAD_X(tunnel_mass_p)
    NITCAD_X(tunnel_mass_evb)
#undef NITCAD_X
    return nullptr;
}

// --- the builder -------------------------------------------------------------------------------

enum class MeshType { tensor, generate2d, generate3d, extrude3d, file };

// lo and hi per dimension, inclusive (DK8).
struct Box {
    std::vector<double> bounds;  // lo0 hi0 lo1 hi1 ...
    bool contains(const mesh::Point& p) const {
        for (std::size_t k = 0; 2 * k < bounds.size(); ++k) {
            if (p[k] < bounds[2 * k] || p[k] > bounds[2 * k + 1]) return false;
        }
        return true;
    }
};

struct RegionSpec {
    const Section* section = nullptr;
    std::optional<device::Region> region;
    std::vector<const Entry*> boxes, polygons, seeds;
    std::vector<const Entry*> extrudes;  // Unit 53 (FD13): extrude3d meshes
};

class Builder {
public:
    explicit Builder(const Deck& d) : d_(d) {}

    std::expected<Plan, base::Error> run();

private:
    base::Error at(int line, const std::string& what) const {
        return {base::ErrorCode::invalid_input, d_.source + ":" + std::to_string(line) + ": " + what,
                base::ErrorContext{.index = static_cast<std::size_t>(line), .value = std::nullopt}};
    }
    // A library error, located: its code and context kept, the message prefixed.
    base::Error located(int line, const std::string& where, base::Error e) const {
        e.message = d_.source + ":" + std::to_string(line) + ": " + where + ": " + e.message;
        return e;
    }
    template <class T>
    std::expected<T, base::Error> lift(int line, std::expected<T, std::string> v) const {
        if (!v) return std::unexpected(at(line, v.error()));
        return std::move(*v);
    }
    std::expected<double, base::Error> num(const Entry& e) const { return lift(e.line, number(e.value)); }
    std::expected<double, base::Error> len(const Entry& e) const { return lift(e.line, length(e.value)); }
    std::expected<int, base::Error> count(const Entry& e, long long lo, long long hi) const {
        auto v = integer(e.value);
        if (!v) return std::unexpected(at(e.line, v.error()));
        if (*v < lo || *v > hi) {
            return std::unexpected(at(e.line, "'" + e.key + "' must lie in " + std::to_string(lo) +
                                                  ".." + std::to_string(hi)));
        }
        return static_cast<int>(*v);
    }
    std::expected<bool, base::Error> flag(const Entry& e) const { return lift(e.line, boolean(e.value)); }
    std::filesystem::path file(const Entry& e) const { return d_.directory / unquote(e.value); }
    base::Error unknown(const Entry& e, const Section& s) const {
        return at(e.line, "unknown key '" + e.key + "' in [" + s.kind + "]");
    }
    // `n` lengths separated by blanks or commas.
    std::expected<std::vector<double>, base::Error> lengths(const Entry& e, std::string_view text,
                                                            std::size_t n) const {
        std::vector<double> out;
        for (const std::string_view item : items(text)) {
            for (const std::string_view t : tokens(item)) {
                NITCAD_DECK_TRY(v, lift(e.line, length(t)));
                out.push_back(v);
            }
        }
        if (out.size() != n) {
            return std::unexpected(at(e.line, "'" + e.key + "' needs " + std::to_string(n) +
                                                  " lengths, has " + std::to_string(out.size())));
        }
        return out;
    }
    // Points of `dim` coordinates separated by commas.
    std::expected<std::vector<mesh::Point>, base::Error> points(const Entry& e, std::string_view text,
                                                                int dim) const {
        std::vector<mesh::Point> out;
        for (const std::string_view item : items(text)) {
            const auto t = tokens(item);
            if (t.size() != static_cast<std::size_t>(dim)) {
                return std::unexpected(at(e.line, "a point needs " + std::to_string(dim) +
                                                      " coordinates: '" + std::string(item) + "'"));
            }
            mesh::Point p{0, 0, 0};
            for (int k = 0; k < dim; ++k) {
                NITCAD_DECK_TRY(v, lift(e.line, length(t[static_cast<std::size_t>(k)])));
                p[static_cast<std::size_t>(k)] = v;
            }
            out.push_back(p);
        }
        return out;
    }

    const Section* single(std::string_view kind) const {
        for (const Section& s : d_.sections) {
            if (s.kind == kind) return &s;
        }
        return nullptr;
    }

    std::expected<void, base::Error> check_sections() const;
    std::expected<RegionSpec, base::Error> region(const Section& s) const;
    std::expected<std::vector<double>, base::Error> axis(const Entry& e) const;
    std::expected<void, base::Error> build_mesh(device::DeviceDescription& d);
    std::expected<void, base::Error> tensor_mesh(const Section& s, device::DeviceDescription& d);
    std::expected<void, base::Error> generated_2d(const Section& s);
    std::expected<void, base::Error> generated_3d(const Section& s);
    std::expected<void, base::Error> extruded_3d(const Section& s);
    // The meshing keys generate3d and extrude3d share; true when `e` was one of them.
    std::expected<bool, base::Error> option_3d(const Entry& e, mesh::MeshOptions3D& o) const;
    std::expected<void, base::Error> file_mesh(const Section& s, device::DeviceDescription& d);
    std::expected<void, base::Error> doping(device::DeviceDescription& d) const;
    std::expected<device::Contact, base::Error> contact(const Section& s,
                                                        const device::DeviceDescription& d) const;
    std::expected<device::Interface, base::Error> interface(const Section& s) const;
    std::expected<device::ThermalContact, base::Error> thermal_contact(const Section& s) const;
    std::expected<void, base::Error> models(assemble::PhysicsModels& m) const;
    std::expected<void, base::Error> solver(Plan& plan) const;
    std::expected<void, base::Error> analysis(Plan& plan) const;
    std::expected<void, base::Error> transient(const Section& s, Plan& plan) const;
    std::expected<void, base::Error> small_signal(const Section& s, Plan& plan) const;
    std::expected<void, base::Error> trace(const Section& s, Plan& plan) const;
    template <class Other>
    std::expected<void, base::Error> operating_points(const Section& s, Plan& plan, Other&& other) const;
    std::expected<bool, base::Error> held_temperature(const Entry& e, Plan& plan,
                                                      std::vector<bool>& given) const;
    std::expected<std::size_t, base::Error> contact_index(const Entry& e, std::string_view name,
                                                          const Plan& plan) const;
    std::expected<std::size_t, base::Error> thermal_index(const Entry& e, std::string_view name,
                                                          const Plan& plan) const;
    base::Error no_restart(const Entry& e, const Section& s) const;
    std::expected<void, base::Error> outputs(Plan& plan) const;
    std::expected<Box, base::Error> box(const Entry& e, int dim) const;

    const Deck& d_;
    std::vector<RegionSpec> regions_;
    MeshType type_ = MeshType::tensor;
    std::optional<mesh::GeneratedMesh> generated_;
};

std::expected<void, base::Error> Builder::check_sections() const {
    static const std::set<std::string_view> once{"device",       "mesh",      "doping",
                                                 "models",       "solver",    "equilibrium",
                                                 "sweep",        "transient", "small_signal",
                                                 "trace",        "output"};
    static const std::set<std::string_view> many{"region", "contact", "interface",
                                                 "thermal_contact"};
    std::set<std::string> seen;
    for (const Section& s : d_.sections) {
        if (once.contains(s.kind)) {
            if (!s.name.empty()) return std::unexpected(at(s.line, "[" + s.kind + "] takes no name"));
            if (!seen.insert(s.kind).second) {
                return std::unexpected(at(s.line, "a second [" + s.kind + "] section"));
            }
        } else if (many.contains(s.kind)) {
            if (s.name.empty()) return std::unexpected(at(s.line, "[" + s.kind + "] needs a name"));
        } else {
            return std::unexpected(at(s.line, "unknown section [" + s.kind + "]"));
        }
    }
    if (!seen.contains("mesh")) return std::unexpected(at(1, "the deck has no [mesh] section"));
    const auto analyses = std::ranges::count_if(
        std::array{"equilibrium", "sweep", "transient", "small_signal", "trace"},
        [&](const char* k) { return seen.contains(k); });
    if (analyses != 1) {
        return std::unexpected(at(1, "the deck needs exactly one analysis: [equilibrium], [sweep], "
                                     "[transient], [small_signal] or [trace]"));
    }
    return {};
}

std::expected<RegionSpec, base::Error> Builder::region(const Section& s) const {
    RegionSpec spec;
    spec.section = &s;
    const Entry* material = nullptr;
    for (const Entry& e : s.entries) {
        if (e.key == "material") {
            if (material) return std::unexpected(at(e.line, "a second 'material'"));
            material = &e;
        }
    }
    if (!material) return std::unexpected(at(s.line, "[region " + s.name + "] has no 'material'"));
    const auto words = tokens(material->value);
    const std::string_view kind = words.front();
    std::variant<physics::SemiconductorParameters, physics::InsulatorParameters> params;
    if (kind == "algaas") {
        if (words.size() < 2 || words.size() > 3) {
            return std::unexpected(at(material->line, "algaas takes the Al fraction and optionally the conduction share"));
        }
        NITCAD_DECK_TRY(x, lift(material->line, number(words[1])));
        double share = 0.85 / 1.247;
        if (words.size() == 3) {
            NITCAD_DECK_TRY(v, lift(material->line, number(words[2])));
            share = v;
        }
        auto p = physics::algaas_parameters(x, share);
        if (!p) return std::unexpected(located(material->line, "material", std::move(p.error())));
        params = *p;
    } else {
        if (words.size() != 1) return std::unexpected(at(material->line, "'" + std::string(kind) + "' takes no arguments"));
        if (kind == "silicon") params = physics::silicon_parameters;
        else if (kind == "silicon_legacy") params = physics::legacy_silicon_parameters;
        else if (kind == "germanium") params = physics::germanium_parameters;
        else if (kind == "gaas") params = physics::gallium_arsenide_parameters;
        else if (kind == "ingaas") params = physics::indium_gallium_arsenide_parameters;
        else if (kind == "sic_4h") params = physics::silicon_carbide_4h_parameters;
        else if (kind == "sio2") params = physics::silicon_dioxide_parameters;
        else if (kind == "insulator") params = physics::InsulatorParameters{.eps_r = 0.0};
        else return std::unexpected(at(material->line, "unknown material '" + std::string(kind) + "'"));
    }
    for (const Entry& e : s.entries) {
        if (&e == material) continue;
        if (e.key == "box") { spec.boxes.push_back(&e); continue; }
        if (e.key == "polygon") { spec.polygons.push_back(&e); continue; }
        if (e.key == "seed") { spec.seeds.push_back(&e); continue; }
        if (e.key == "extrude") { spec.extrudes.push_back(&e); continue; }
        if (auto* sp = std::get_if<physics::SemiconductorParameters>(&params)) {
            if (e.key == "band_to_band.direct_gap") {
                NITCAD_DECK_TRY(v, flag(e));
                sp->band_to_band.direct_gap = v;
                continue;
            }
            if (double* target = semiconductor_parameter(*sp, e.key)) {
                NITCAD_DECK_TRY(v, num(e));
                *target = v;
                continue;
            }
        } else if (double* target = insulator_parameter(std::get<physics::InsulatorParameters>(params), e.key)) {
            NITCAD_DECK_TRY(v, num(e));
            *target = v;
            continue;
        }
        return std::unexpected(at(e.line, "unknown key '" + e.key + "' in [region " + s.name +
                                              "] (not a parameter of its material)"));
    }
    if (const auto* sp = std::get_if<physics::SemiconductorParameters>(&params)) {
        auto m = physics::Semiconductor::create(*sp);
        if (!m) return std::unexpected(located(material->line, "[region " + s.name + "]", std::move(m.error())));
        spec.region = {s.name, std::move(*m)};
    } else {
        auto m = physics::Insulator::create(std::get<physics::InsulatorParameters>(params));
        if (!m) return std::unexpected(located(material->line, "[region " + s.name + "]", std::move(m.error())));
        spec.region = {s.name, std::move(*m)};
    }
    return spec;
}

// Pieces joined by '+': linspace(a, b, n) or a comma list of lengths. A '+' right after an
// exponent's 'e' belongs to the number.
std::expected<std::vector<double>, base::Error> Builder::axis(const Entry& e) const {
    std::vector<std::string_view> pieces;
    const std::string_view v = e.value;
    std::size_t start = 0;
    int depth = 0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (v[i] == '(') ++depth;
        if (v[i] == ')') --depth;
        const bool exponent = i > 0 && (v[i - 1] == 'e' || v[i - 1] == 'E');
        if (v[i] == '+' && depth == 0 && !exponent) {
            pieces.push_back(trimmed(v.substr(start, i - start)));
            start = i + 1;
        }
    }
    pieces.push_back(trimmed(v.substr(start)));
    std::vector<double> out;
    const auto append = [&out](double x) {
        if (out.empty() || out.back() != x) out.push_back(x);
    };
    for (const std::string_view piece : pieces) {
        if (piece.starts_with("linspace(")) {
            if (!piece.ends_with(")")) return std::unexpected(at(e.line, "linspace( without ')'"));
            const auto args = items(piece.substr(9, piece.size() - 10));
            if (args.size() != 3) return std::unexpected(at(e.line, "linspace takes (first, last, count)"));
            NITCAD_DECK_TRY(a, lift(e.line, length(args[0])));
            NITCAD_DECK_TRY(b, lift(e.line, length(args[1])));
            NITCAD_DECK_TRY(n, lift(e.line, integer(args[2])));
            if (n < 2 || n > 10'000'000) return std::unexpected(at(e.line, "linspace needs a count of at least 2"));
            for (long long k = 0; k < n; ++k) {
                append(k == n - 1 ? b : a + (b - a) * static_cast<double>(k) / static_cast<double>(n - 1));
            }
            continue;
        }
        for (const std::string_view item : items(piece)) {
            NITCAD_DECK_TRY(x, lift(e.line, length(item)));
            append(x);
        }
    }
    return out;
}

std::expected<Box, base::Error> Builder::box(const Entry& e, int dim) const {
    NITCAD_DECK_TRY(b, lengths(e, e.value, static_cast<std::size_t>(2 * dim)));
    for (int k = 0; k < dim; ++k) {
        if (b[static_cast<std::size_t>(2 * k)] > b[static_cast<std::size_t>(2 * k + 1)]) {
            return std::unexpected(at(e.line, "a box needs lo <= hi on each axis"));
        }
    }
    return Box{std::move(b)};
}

std::expected<void, base::Error> Builder::tensor_mesh(const Section& s, device::DeviceDescription& d) {
    std::array<std::optional<std::vector<double>>, 3> axes;
    std::vector<const Entry*> straddles;
    static constexpr std::array<std::string_view, 3> names{"x", "y", "z"};
    for (const Entry& e : s.entries) {
        if (e.key == "type") continue;
        if (e.key.starts_with("straddle ")) {
            straddles.push_back(&e);
            continue;
        }
        const auto it = std::find(names.begin(), names.end(), e.key);
        if (it == names.end()) return std::unexpected(unknown(e, s));
        auto& a = axes[static_cast<std::size_t>(it - names.begin())];
        if (a) return std::unexpected(at(e.line, "a second '" + e.key + "'"));
        NITCAD_DECK_TRY(v, axis(e));
        a = std::move(v);
    }
    if (!axes[0] || (axes[2] && !axes[1])) {
        return std::unexpected(at(s.line, "a tensor mesh needs x, then y, then z"));
    }
    const int dim = axes[2] ? 3 : axes[1] ? 2 : 1;
    for (const Entry* e : straddles) {
        const auto which = e->key.substr(9);
        const auto it = std::find(names.begin(), names.begin() + dim, which);
        if (it == names.begin() + dim) return std::unexpected(at(e->line, "straddle on an axis the mesh does not have"));
        NITCAD_DECK_TRY(ps, lengths(*e, e->value, 2));
        auto& a = *axes[static_cast<std::size_t>(it - names.begin())];
        auto moved = mesh::straddle_interface(std::move(a), ps[0], ps[1]);
        if (!moved) return std::unexpected(located(e->line, "straddle", std::move(moved.error())));
        a = std::move(*moved);
    }
    auto m = dim == 1   ? mesh::make_tensor_grid(*axes[0])
             : dim == 2 ? mesh::make_tensor_grid(*axes[0], *axes[1])
                        : mesh::make_tensor_grid(*axes[0], *axes[1], *axes[2]);
    if (!m) return std::unexpected(located(s.line, "[mesh]", std::move(m.error())));
    auto cells = dim == 1   ? mesh::TensorCells::create(*m, *axes[0])
                 : dim == 2 ? mesh::TensorCells::create(*m, *axes[0], *axes[1])
                            : mesh::TensorCells::create(*m, *axes[0], *axes[1], *axes[2]);
    if (!cells) return std::unexpected(located(s.line, "[mesh]", std::move(cells.error())));
    // Regions: the background (the one region without boxes) first, then every box in deck order.
    std::optional<device::RegionId> background;
    for (std::size_t r = 0; r < regions_.size(); ++r) {
        const RegionSpec& g = regions_[r];
        if (!g.polygons.empty() || !g.seeds.empty()) {
            return std::unexpected(at((g.polygons.empty() ? g.seeds : g.polygons).front()->line,
                                      "polygon and seed are for generated meshes; a tensor mesh takes box"));
        }
        if (g.boxes.empty()) {
            if (background) return std::unexpected(at(g.section->line, "two regions without a box: only one may be the background"));
            background = static_cast<device::RegionId>(r);
        }
    }
    constexpr device::RegionId none = -1;
    d.node_region.assign(m->node_count(), background.value_or(none));
    for (std::size_t r = 0; r < regions_.size(); ++r) {
        for (const Entry* e : regions_[r].boxes) {
            NITCAD_DECK_TRY(b, box(*e, dim));
            for (std::size_t i = 0; i < m->node_count(); ++i) {
                if (b.contains(m->points()[i])) d.node_region[i] = static_cast<device::RegionId>(r);
            }
        }
    }
    const auto unassigned = std::find(d.node_region.begin(), d.node_region.end(), none);
    if (unassigned != d.node_region.end()) {
        const auto& p = m->points()[static_cast<std::size_t>(unassigned - d.node_region.begin())];
        return std::unexpected(at(s.line, "the node at (" + std::to_string(p[0]) + ", " + std::to_string(p[1]) +
                                              ", " + std::to_string(p[2]) + ") cm is in no region's box"));
    }
    d.mesh = std::move(*m);
    d.cells = std::move(*cells);
    return {};
}

std::expected<void, base::Error> Builder::generated_2d(const Section& s) {
    mesh::Geometry2D g;
    mesh::MeshOptions2D o;
    for (const Entry& e : s.entries) {
        if (e.key == "type") continue;
        if (e.key.starts_with("patch ")) {
            NITCAD_DECK_TRY(ab, points(e, e.value, 2));
            if (ab.size() != 2) return std::unexpected(at(e.line, "a patch is two points"));
            g.patches.push_back({e.key.substr(6), ab[0], ab[1]});
        } else if (e.key == "max_size") {
            NITCAD_DECK_TRY(v, len(e));
            o.max_size = v;
        } else if (e.key == "interface_resolution") {
            NITCAD_DECK_TRY(v, len(e));
            o.interface_resolution = v;
        } else if (e.key == "min_angle") {
            NITCAD_DECK_TRY(v, num(e));
            o.min_angle_deg = v;
        } else if (e.key == "max_nodes") {
            NITCAD_DECK_TRY(v, count(e, 1, std::numeric_limits<int>::max()));
            o.max_nodes = static_cast<std::size_t>(v);
        } else if (e.key == "box") {
            NITCAD_DECK_TRY(b, lengths(e, e.value, 5));
            o.boxes.push_back({{b[0], b[1], 0.0}, {b[2], b[3], 0.0}, b[4]});
        } else {
            return std::unexpected(unknown(e, s));
        }
    }
    for (std::size_t r = 0; r < regions_.size(); ++r) {
        const RegionSpec& spec = regions_[r];
        if (!spec.boxes.empty() || !spec.seeds.empty()) {
            return std::unexpected(at((spec.boxes.empty() ? spec.seeds : spec.boxes).front()->line,
                                      "a generate2d mesh takes polygon in its regions"));
        }
        for (const Entry* e : spec.polygons) {
            NITCAD_DECK_TRY(v, points(*e, e->value, 2));
            g.polygons.push_back({std::move(v), static_cast<std::int32_t>(r)});
        }
    }
    auto m = mesh::generate_mesh(g, o);
    if (!m) return std::unexpected(located(s.line, "[mesh] generate2d", std::move(m.error())));
    generated_ = std::move(*m);
    return {};
}

std::expected<void, base::Error> Builder::generated_3d(const Section& s) {
    mesh::Geometry3D g;
    mesh::MeshOptions3D o;
    std::vector<const Entry*> patches;
    bool cuboid = false;
    for (const Entry& e : s.entries) {
        if (e.key == "type") continue;
        if (e.key == "vertex") {
            NITCAD_DECK_TRY(p, points(e, e.value, 3));
            if (p.size() != 1) return std::unexpected(at(e.line, "a vertex is one point"));
            g.domain.vertices.push_back(p[0]);
        } else if (e.key.starts_with("face ")) {
            mesh::PolyhedronFace f;
            f.name = e.key.substr(5);
            for (const std::string_view t : tokens(e.value)) {
                NITCAD_DECK_TRY(k, lift(e.line, integer(t)));
                if (k < 0) return std::unexpected(at(e.line, "a vertex index is >= 0"));
                f.vertices.push_back(static_cast<std::size_t>(k));
            }
            g.domain.faces.push_back(std::move(f));
        } else if (e.key == "cuboid") {
            NITCAD_DECK_TRY(b, lengths(e, e.value, 6));
            cuboid = true;
            for (int k = 0; k < 8; ++k) {
                g.domain.vertices.push_back({b[(k & 1) ? 1 : 0], b[(k & 2) ? 3 : 2], b[(k & 4) ? 5 : 4]});
            }
            // Vertex k has bit a set for the upper end along axis a; each face listed around itself.
            g.domain.faces = {{{0, 4, 6, 2}, "x_min", {}}, {{1, 3, 7, 5}, "x_max", {}},
                              {{0, 1, 5, 4}, "y_min", {}}, {{2, 6, 7, 3}, "y_max", {}},
                              {{0, 2, 3, 1}, "z_min", {}}, {{4, 5, 7, 6}, "z_max", {}}};
        } else if (e.key.starts_with("patch ")) {
            patches.push_back(&e);
        } else if (e.key == "plane") {
            NITCAD_DECK_TRY(p, points(e, e.value, 3));
            if (p.size() != 3) return std::unexpected(at(e.line, "a plane is three points"));
            g.interfaces.push_back({{p[0], p[1], p[2]}});
        } else {
            NITCAD_DECK_TRY(taken, option_3d(e, o));
            if (!taken) return std::unexpected(unknown(e, s));
        }
    }
    if (cuboid && g.domain.faces.size() != 6) {
        return std::unexpected(at(s.line, "cuboid replaces vertex and face: use one or the other"));
    }
    // patch <face> <name> = points: a tiling of the face into named patches.
    for (const Entry* e : patches) {
        const auto w = tokens(e->key);
        if (w.size() != 3) return std::unexpected(at(e->line, "a face patch is 'patch <face> <name>'"));
        auto f = std::find_if(g.domain.faces.begin(), g.domain.faces.end(),
                              [&](const mesh::PolyhedronFace& x) { return x.name == w[1]; });
        if (f == g.domain.faces.end()) return std::unexpected(at(e->line, "no face named '" + std::string(w[1]) + "'"));
        NITCAD_DECK_TRY(p, points(*e, e->value, 3));
        f->patches.push_back({std::string(w[2]), std::move(p)});
    }
    for (std::size_t r = 0; r < regions_.size(); ++r) {
        const RegionSpec& spec = regions_[r];
        if (!spec.boxes.empty() || !spec.polygons.empty()) {
            return std::unexpected(at((spec.boxes.empty() ? spec.polygons : spec.boxes).front()->line,
                                      "a generate3d mesh takes seed in its regions"));
        }
        for (const Entry* e : spec.seeds) {
            NITCAD_DECK_TRY(p, points(*e, e->value, 3));
            if (p.size() != 1) return std::unexpected(at(e->line, "a seed is one point"));
            g.regions.push_back({p[0], static_cast<std::int32_t>(r)});
        }
    }
    if (g.interfaces.empty()) {
        if (regions_.size() != 1) return std::unexpected(at(s.line, "a generate3d mesh without planes has one region"));
        g.region = 0;
        g.regions.clear();
    }
    auto m = mesh::generate_mesh(g, o);
    if (!m) return std::unexpected(located(s.line, "[mesh] generate3d", std::move(m.error())));
    generated_ = std::move(*m);
    return {};
}

std::expected<bool, base::Error> Builder::option_3d(const Entry& e, mesh::MeshOptions3D& o) const {
    if (e.key == "max_size") {
        NITCAD_DECK_TRY(v, len(e));
        o.max_size = v;
    } else if (e.key == "interface_resolution") {
        NITCAD_DECK_TRY(v, len(e));
        o.interface_resolution = v;
    } else if (e.key == "radius_edge") {
        NITCAD_DECK_TRY(v, num(e));
        o.radius_edge = v;
    } else if (e.key == "max_nodes") {
        NITCAD_DECK_TRY(v, count(e, 1, std::numeric_limits<int>::max()));
        o.max_nodes = static_cast<std::size_t>(v);
    } else if (e.key == "box") {
        NITCAD_DECK_TRY(b, lengths(e, e.value, 7));
        o.boxes.push_back({{b[0], b[1], b[2]}, {b[3], b[4], b[5]}, b[6]});
    } else {
        return false;
    }
    return true;
}

// Unit 53 (DECISIONS.md FD13): a sectioned extrusion. `[mesh]` gives the axis, the stations, the
// cells and the patch boxes; each region assigns its cells per section with `extrude`.
std::expected<void, base::Error> Builder::extruded_3d(const Section& s) {
    mesh::ExtrudedGeometry x;
    mesh::MeshOptions3D o;
    const Entry* axis_entry = nullptr;
    const Entry* stations_entry = nullptr;
    for (const Entry& e : s.entries) {
        if (e.key == "type") continue;
        if (e.key == "axis") {
            if (e.value != "x" && e.value != "y" && e.value != "z") return std::unexpected(at(e.line, "axis is x, y or z"));
            axis_entry = &e;
            x.axis = e.value == "x" ? 0 : e.value == "y" ? 1 : 2;
        } else if (e.key == "stations") {
            stations_entry = &e;
            for (const std::string_view item : items(e.value)) {
                for (const std::string_view t : tokens(item)) {
                    NITCAD_DECK_TRY(v, lift(e.line, length(t)));
                    x.stations.push_back(v);
                }
            }
        } else if (e.key.starts_with("cell ")) {
            const auto w = tokens(e.key);
            if (w.size() != 2) return std::unexpected(at(e.line, "a cell is 'cell <name>'"));
            NITCAD_DECK_TRY(polygon, points(e, e.value, 2));
            x.cells.push_back({std::string(w[1]), std::move(polygon)});
        } else if (e.key.starts_with("patch ")) {
            const auto w = tokens(e.key);
            if (w.size() != 2) return std::unexpected(at(e.line, "a patch is 'patch <name> = <x0 y0 z0 x1 y1 z1>'"));
            NITCAD_DECK_TRY(b, lengths(e, e.value, 6));
            if (b[0] > b[3] || b[1] > b[4] || b[2] > b[5]) return std::unexpected(at(e.line, "a box needs lo <= hi on each axis"));
            x.patches.push_back({std::string(w[1]), {b[0], b[1], b[2]}, {b[3], b[4], b[5]}});
        } else {
            NITCAD_DECK_TRY(taken, option_3d(e, o));
            if (!taken) return std::unexpected(unknown(e, s));
        }
    }
    if (!axis_entry) return std::unexpected(at(s.line, "an extrude3d mesh needs 'axis'"));
    if (!stations_entry) return std::unexpected(at(s.line, "an extrude3d mesh needs 'stations'"));
    if (x.cells.empty()) return std::unexpected(at(s.line, "an extrude3d mesh needs at least one 'cell'"));
    const std::size_t sections = x.stations.size() > 1 ? x.stations.size() - 1 : 0;
    x.regions.assign(x.cells.size(), std::vector<std::int32_t>(sections, mesh::ExtrudedGeometry::void_region));
    // region extrude = <cells> <t_lo>:<t_hi>; the range runs between two stations.
    for (std::size_t r = 0; r < regions_.size(); ++r) {
        const RegionSpec& spec = regions_[r];
        if (!spec.boxes.empty() || !spec.polygons.empty() || !spec.seeds.empty()) {
            return std::unexpected(at((!spec.boxes.empty() ? spec.boxes : !spec.polygons.empty() ? spec.polygons : spec.seeds).front()->line,
                                      "an extrude3d mesh takes extrude in its regions"));
        }
        for (const Entry* e : spec.extrudes) {
            const auto w = tokens(e->value);
            if (w.size() < 2) return std::unexpected(at(e->line, "extrude is '<cells> <t_lo>:<t_hi>'"));
            const auto range = items(w.back(), ':');
            if (range.size() != 2) return std::unexpected(at(e->line, "the range is <t_lo>:<t_hi>"));
            NITCAD_DECK_TRY(lo, lift(e->line, length(range[0])));
            NITCAD_DECK_TRY(hi, lift(e->line, length(range[1])));
            const auto station = [&](double t) -> std::optional<std::size_t> {
                const auto it = std::find(x.stations.begin(), x.stations.end(), t);
                if (it == x.stations.end()) return std::nullopt;
                return static_cast<std::size_t>(it - x.stations.begin());
            };
            const auto a = station(lo);
            const auto b = station(hi);
            if (!a || !b) return std::unexpected(at(e->line, "the range's ends must be stations"));
            if (*a >= *b) return std::unexpected(at(e->line, "the range needs t_lo < t_hi"));
            for (std::size_t k = 0; k + 1 < w.size(); ++k) {
                const auto c = std::find_if(x.cells.begin(), x.cells.end(), [&](const mesh::ExtrusionCell& cell) { return cell.name == w[k]; });
                if (c == x.cells.end()) return std::unexpected(at(e->line, "no cell named '" + std::string(w[k]) + "'"));
                auto& row = x.regions[static_cast<std::size_t>(c - x.cells.begin())];
                for (std::size_t j = *a; j < *b; ++j) {
                    if (row[j] != mesh::ExtrudedGeometry::void_region) {
                        return std::unexpected(at(e->line, "cell '" + std::string(w[k]) + "' is already assigned in section " + std::to_string(j)));
                    }
                    row[j] = static_cast<std::int32_t>(r);
                }
            }
        }
    }
    auto complex = mesh::extrude(x);
    if (!complex) return std::unexpected(located(s.line, "[mesh] extrude3d", std::move(complex.error())));
    auto m = mesh::generate_mesh(*complex, o);
    if (!m) return std::unexpected(located(s.line, "[mesh] extrude3d", std::move(m.error())));
    generated_ = std::move(*m);
    return {};
}

std::expected<void, base::Error> Builder::file_mesh(const Section& s, device::DeviceDescription& d) {
    const Entry* path = nullptr;
    for (const Entry& e : s.entries) {
        if (e.key == "type") continue;
        if (e.key != "path") return std::unexpected(unknown(e, s));
        path = &e;
    }
    if (!path) return std::unexpected(at(s.line, "a file mesh needs 'path'"));
    for (const RegionSpec& g : regions_) {
        if (!g.boxes.empty() || !g.polygons.empty() || !g.seeds.empty()) {
            return std::unexpected(at(g.section->line, "a file mesh carries its regions: no box, polygon or seed"));
        }
    }
    auto generated = io::load_generated_mesh(file(*path));
    if (generated) {
        generated_ = std::move(*generated);
        return {};
    }
    if (generated.error().message.find("plain mesh") == std::string::npos) {
        return std::unexpected(located(path->line, "[mesh] path", std::move(generated.error())));
    }
    auto m = io::load_mesh(file(*path));
    if (!m) return std::unexpected(located(path->line, "[mesh] path", std::move(m.error())));
    if (regions_.size() != 1) return std::unexpected(at(path->line, "a plain mesh file has one region"));
    d.node_region.assign(m->node_count(), 0);
    d.mesh = std::move(*m);
    return {};
}

std::expected<void, base::Error> Builder::build_mesh(device::DeviceDescription& d) {
    const Section& s = *single("mesh");
    const Entry* type = nullptr;
    for (const Entry& e : s.entries) {
        if (e.key == "type") type = &e;
    }
    if (!type) return std::unexpected(at(s.line, "[mesh] needs 'type'"));
    if (type->value == "tensor") type_ = MeshType::tensor;
    else if (type->value == "generate2d") type_ = MeshType::generate2d;
    else if (type->value == "generate3d") type_ = MeshType::generate3d;
    else if (type->value == "extrude3d") type_ = MeshType::extrude3d;
    else if (type->value == "file") type_ = MeshType::file;
    else return std::unexpected(at(type->line, "unknown mesh type '" + type->value + "'"));
    if (type_ != MeshType::extrude3d) {
        for (const RegionSpec& g : regions_) {
            if (!g.extrudes.empty()) return std::unexpected(at(g.extrudes.front()->line, "extrude is for an extrude3d mesh"));
        }
    }
    std::expected<void, base::Error> made;
    switch (type_) {
    case MeshType::extrude3d: made = extruded_3d(s); break;
    case MeshType::tensor: made = tensor_mesh(s, d); break;
    case MeshType::generate2d: made = generated_2d(s); break;
    case MeshType::generate3d: made = generated_3d(s); break;
    case MeshType::file: made = file_mesh(s, d); break;
    }
    if (!made) return made;
    if (generated_) {
        d.mesh = generated_->mesh;
        d.node_region.assign(generated_->node_region.begin(), generated_->node_region.end());
        d.interface_certificate = generated_->certificate;
        for (const device::RegionId r : d.node_region) {
            if (r < 0 || static_cast<std::size_t>(r) >= regions_.size()) {
                return std::unexpected(at(s.line, "the mesh has region " + std::to_string(r) + " but the deck declares " +
                                                      std::to_string(regions_.size()) + " region(s)"));
            }
        }
    }
    return {};
}

std::expected<void, base::Error> Builder::doping(device::DeviceDescription& d) const {
    std::vector<device::DopingProfile> profiles;
    if (const Section* s = single("doping")) {
        for (const Entry& e : s->entries) {
            if (e.key != "profile") return std::unexpected(unknown(e, *s));
            const auto t = tokens(e.value);
            if (t.size() < 3) return std::unexpected(at(e.line, "a profile starts '<donor|acceptor> <shape> <peak>'"));
            device::DopingProfile p;
            if (t[0] == "donor") p.dopant = device::Dopant::donor;
            else if (t[0] == "acceptor") p.dopant = device::Dopant::acceptor;
            else return std::unexpected(at(e.line, "the dopant is donor or acceptor"));
            if (t[1] == "constant") p.shape = device::DepthShape::constant;
            else if (t[1] == "gaussian") p.shape = device::DepthShape::gaussian;
            else if (t[1] == "erfc") p.shape = device::DepthShape::erfc;
            else return std::unexpected(at(e.line, "the shape is constant, gaussian or erfc"));
            NITCAD_DECK_TRY(peak, lift(e.line, number(t[2])));
            p.peak_cm3 = peak;
            for (std::size_t k = 3; k < t.size(); ++k) {
                const std::size_t eq = t[k].find('=');
                if (eq == std::string_view::npos) return std::unexpected(at(e.line, "expected name=value, not '" + std::string(t[k]) + "'"));
                const std::string_view name = t[k].substr(0, eq);
                const std::string_view value = t[k].substr(eq + 1);
                static constexpr std::array<std::string_view, 3> axes{"x", "y", "z"};
                const auto a = std::find(axes.begin(), axes.end(), name);
                if (name == "axis") {
                    const auto v = std::find(axes.begin(), axes.end(), value);
                    if (v == axes.end()) return std::unexpected(at(e.line, "axis is x, y or z"));
                    p.axis = static_cast<int>(v - axes.begin());
                } else if (name == "center" || name == "width" || name == "lateral") {
                    NITCAD_DECK_TRY(v, lift(e.line, length(value)));
                    (name == "center" ? p.center : name == "width" ? p.width : p.lateral) = v;
                } else if (a != axes.end()) {
                    const auto range = items(value, ':');
                    if (range.size() != 2) return std::unexpected(at(e.line, "a window is <lo>:<hi>"));
                    NITCAD_DECK_TRY(lo, lift(e.line, length(range[0])));
                    NITCAD_DECK_TRY(hi, lift(e.line, length(range[1])));
                    const auto k3 = static_cast<std::size_t>(a - axes.begin());
                    p.lo[k3] = lo;
                    p.hi[k3] = hi;
                } else {
                    return std::unexpected(at(e.line, "unknown profile attribute '" + std::string(name) + "'"));
                }
            }
            profiles.push_back(p);
        }
    }
    auto sampled = device::sample(profiles, d.mesh, d.node_region, d.regions);
    if (!sampled) {
        const Section* s = single("doping");
        const int line = s && sampled.error().context && sampled.error().context->index &&
                                 *sampled.error().context->index < s->entries.size()
                             ? s->entries[*sampled.error().context->index].line
                             : (s ? s->line : 1);
        return std::unexpected(located(line, "[doping]", std::move(sampled.error())));
    }
    d.donors = std::move(sampled->donors);
    d.acceptors = std::move(sampled->acceptors);
    return {};
}

std::expected<device::Contact, base::Error> Builder::contact(const Section& s,
                                                             const device::DeviceDescription& d) const {
    device::Contact c{.name = s.name, .kind = device::ContactKind::ohmic, .nodes = {}};
    const Entry* patch = nullptr;
    const Entry* kind = nullptr;
    std::vector<const Entry*> boxes;
    bool gate_keys = false;
    const Entry* work_function = nullptr;
    const Entry* barrier = nullptr;
    bool schottky_keys = false;
    for (const Entry& e : s.entries) {
        if (e.key == "kind") kind = &e;
        else if (e.key == "patch") patch = &e;
        else if (e.key == "box") boxes.push_back(&e);
        else if (e.key == "oxide_thickness") {
            NITCAD_DECK_TRY(v, len(e));
            c.gate.oxide_thickness_cm = v;
            gate_keys = true;
        } else if (e.key == "oxide_permittivity") {
            NITCAD_DECK_TRY(v, num(e));
            c.gate.oxide_relative_permittivity = v;
            gate_keys = true;
        } else if (e.key == "fixed_charge") {
            NITCAD_DECK_TRY(v, num(e));
            c.gate.fixed_charge_cm2 = v;
            gate_keys = true;
        } else if (e.key == "electrode") {
            device::GateElectrode g{};
            if (e.value == "n_poly") g = device::GateElectrode::n_poly;
            else if (e.value == "p_poly") g = device::GateElectrode::p_poly;
            else if (e.value == "metal") g = device::GateElectrode::metal;
            else return std::unexpected(at(e.line, "electrode is n_poly, p_poly or metal"));
            c.gate.electrode = g;
            c.electrode.kind = g;
        } else if (e.key == "work_function") {
            NITCAD_DECK_TRY(v, num(e));
            c.gate.work_function_eV = v;
            c.electrode.work_function_eV = v;
            c.schottky.work_function_eV = v;
            work_function = &e;
        } else if (e.key == "barrier") {  // Unit 48
            NITCAD_DECK_TRY(v, num(e));
            c.schottky.barrier_eV = v;
            barrier = &e;
        } else if (e.key == "richardson_n") {
            NITCAD_DECK_TRY(v, num(e));
            c.schottky.richardson_n = v;
            schottky_keys = true;
        } else if (e.key == "richardson_p") {
            NITCAD_DECK_TRY(v, num(e));
            c.schottky.richardson_p = v;
            schottky_keys = true;
        } else if (e.key == "tunnelling") {  // Unit 51 (GD8): off, on, or processes
            std::uint8_t processes = 0;
            std::string word;
            const auto take = [&]() -> bool {
                if (word.empty()) return true;
                if (word == "off") {
                } else if (word == "on") {
                    processes |= device::all_tunnel_processes;
                } else if (word == "ecb") {
                    processes |= static_cast<std::uint8_t>(device::TunnelProcess::ecb);
                } else if (word == "hvb") {
                    processes |= static_cast<std::uint8_t>(device::TunnelProcess::hvb);
                } else if (word == "evb") {
                    processes |= static_cast<std::uint8_t>(device::TunnelProcess::evb);
                } else {
                    return false;
                }
                word.clear();
                return true;
            };
            for (const char ch : e.value + ",") {
                if (ch == ',' || ch == ' ' || ch == '\t') {
                    if (!take()) {
                        return std::unexpected(at(e.line, "tunnelling is off, on, or a list of ecb, hvb and evb"));
                    }
                } else {
                    word.push_back(ch);
                }
            }
            c.tunnelling = processes;
        } else if (e.key == "oxide_affinity" || e.key == "oxide_band_gap" ||
                   e.key == "oxide_mass_n" || e.key == "oxide_mass_p" ||
                   e.key == "oxide_mass_evb") {  // Unit 51 (GD5)
            NITCAD_DECK_TRY(v, num(e));
            double& target = e.key == "oxide_affinity"   ? c.gate.oxide_affinity_eV
                             : e.key == "oxide_band_gap" ? c.gate.oxide_band_gap_eV
                             : e.key == "oxide_mass_n"   ? c.gate.oxide_mass_n
                             : e.key == "oxide_mass_p"   ? c.gate.oxide_mass_p
                                                         : c.gate.oxide_mass_evb;
            target = v;
            gate_keys = true;
        } else {
            return std::unexpected(unknown(e, s));
        }
    }
    if (!kind) return std::unexpected(at(s.line, "[contact " + s.name + "] needs 'kind'"));
    if (kind->value == "ohmic") c.kind = device::ContactKind::ohmic;
    else if (kind->value == "gate") c.kind = device::ContactKind::gate;
    else if (kind->value == "electrode") c.kind = device::ContactKind::electrode;
    else if (kind->value == "schottky") c.kind = device::ContactKind::schottky;
    else return std::unexpected(at(kind->line, "kind is ohmic, gate, electrode or schottky"));
    if (gate_keys && c.kind != device::ContactKind::gate) {
        return std::unexpected(at(s.line, "oxide_thickness, oxide_permittivity and fixed_charge are for a gate"));
    }
    if ((barrier || schottky_keys) && c.kind != device::ContactKind::schottky) {
        return std::unexpected(at(s.line, "barrier, richardson_n and richardson_p are for a Schottky contact"));
    }
    if (c.kind == device::ContactKind::schottky) {
        // Exactly one of the barrier and the work function (SD1).
        if ((barrier != nullptr) == (work_function != nullptr)) {
            return std::unexpected(at(s.line, "[contact " + s.name + "] needs exactly one of 'barrier' and 'work_function'"));
        }
        c.schottky.by_work_function = work_function != nullptr;
        c.schottky.boundary = patch ? patch->value : std::string{};
    }
    if (!patch) return std::unexpected(at(s.line, "[contact " + s.name + "] needs 'patch'"));
    const mesh::BoundaryPatch* b = d.mesh.find_boundary(patch->value);
    if (!b) return std::unexpected(at(patch->line, "the mesh has no boundary patch '" + patch->value + "'"));
    if (c.kind == device::ContactKind::gate) c.gate.boundary = patch->value;
    const bool want_insulator = c.kind == device::ContactKind::electrode;
    std::vector<Box> parsed;
    for (const Entry* e : boxes) {
        NITCAD_DECK_TRY(x, box(*e, d.mesh.dimension()));
        parsed.push_back(std::move(x));
    }
    for (const mesh::NodeId k : b->nodes) {
        const auto& region = d.regions[static_cast<std::size_t>(d.node_region[static_cast<std::size_t>(k)])];
        if (device::is_insulator(region) != want_insulator) continue;
        const mesh::Point& p = d.mesh.points()[static_cast<std::size_t>(k)];
        // Without boxes the whole patch; with boxes, the nodes inside any of them.
        if (parsed.empty() || std::ranges::any_of(parsed, [&](const Box& x) { return x.contains(p); })) {
            c.nodes.push_back(k);
        }
    }
    std::sort(c.nodes.begin(), c.nodes.end());
    if (c.nodes.empty()) {
        return std::unexpected(at(patch->line, std::string("patch '") + patch->value + "' has no " +
                                                   (want_insulator ? "insulator" : "semiconductor") +
                                                   " node for this contact"));
    }
    return c;
}

std::expected<device::Interface, base::Error> Builder::interface(const Section& s) const {
    const auto names = tokens(s.name);
    if (names.size() != 2) return std::unexpected(at(s.line, "[interface] names two regions"));
    device::Interface f{.region_a = std::string(names[0]), .region_b = std::string(names[1])};
    const auto trap_type = [&](const Entry& e, std::string_view t) -> std::expected<physics::TrapType, base::Error> {
        if (t == "donor") return physics::TrapType::donor;
        if (t == "acceptor") return physics::TrapType::acceptor;
        return std::unexpected(at(e.line, "a trap is donor or acceptor"));
    };
    for (const Entry& e : s.entries) {
        if (e.key == "transport") {
            if (e.value == "drift_diffusion") f.transport = device::InterfaceTransport::drift_diffusion;
            else if (e.value == "thermionic") f.transport = device::InterfaceTransport::thermionic_emission;
            else return std::unexpected(at(e.line, "transport is drift_diffusion or thermionic"));
        } else if (e.key == "fixed_charge") {
            NITCAD_DECK_TRY(v, num(e));
            f.fixed_charge_cm2 = v;
        } else if (e.key == "recombination_velocity_n") {
            NITCAD_DECK_TRY(v, num(e));
            f.recombination_velocity_n_cm_s = v;
        } else if (e.key == "recombination_velocity_p") {
            NITCAD_DECK_TRY(v, num(e));
            f.recombination_velocity_p_cm_s = v;
        } else if (e.key == "thermal_velocity_n") {
            NITCAD_DECK_TRY(v, num(e));
            f.traps.thermal_velocity_n_cm_s = v;
        } else if (e.key == "thermal_velocity_p") {
            NITCAD_DECK_TRY(v, num(e));
            f.traps.thermal_velocity_p_cm_s = v;
        } else if (e.key == "trap_level" || e.key == "trap_band") {
            const auto t = tokens(e.value);
            const bool band = e.key == "trap_band";
            const std::size_t fixed = band ? 4 : 3;
            if (t.size() != fixed && t.size() != fixed + 2) {
                return std::unexpected(at(e.line, band ? "trap_band = <type> <D_it> <E_low> <E_high> [<sigma_n> <sigma_p>]"
                                                       : "trap_level = <type> <N_t> <E> [<sigma_n> <sigma_p>]"));
            }
            NITCAD_DECK_TRY(type, trap_type(e, t[0]));
            std::vector<double> v;
            for (std::size_t k = 1; k < t.size(); ++k) {
                NITCAD_DECK_TRY(x, lift(e.line, number(t[k])));
                v.push_back(x);
            }
            const double sn = t.size() > fixed ? v[fixed - 1] : 1e-15;
            const double sp = t.size() > fixed ? v[fixed] : 1e-15;
            if (band) f.traps.bands.push_back({type, v[0], v[1], v[2], sn, sp});
            else f.traps.levels.push_back({type, v[0], v[1], sn, sp});
        } else {
            return std::unexpected(unknown(e, s));
        }
    }
    return f;
}

std::expected<device::ThermalContact, base::Error> Builder::thermal_contact(const Section& s) const {
    device::ThermalContact t{.name = s.name, .boundary = {}};
    for (const Entry& e : s.entries) {
        if (e.key == "patch") t.boundary = e.value;
        else if (e.key == "kind") {
            if (e.value == "isothermal") t.kind = device::ThermalContactKind::isothermal;
            else if (e.value == "resistance") t.kind = device::ThermalContactKind::resistance;
            else return std::unexpected(at(e.line, "kind is isothermal or resistance"));
        } else if (e.key == "temperature") {
            NITCAD_DECK_TRY(v, num(e));
            t.temperature_K = v;
        } else if (e.key == "resistance") {
            NITCAD_DECK_TRY(v, num(e));
            t.resistance_K_cm2_W = v;
        } else {
            return std::unexpected(unknown(e, s));
        }
    }
    if (t.boundary.empty()) return std::unexpected(at(s.line, "[thermal_contact " + s.name + "] needs 'patch'"));
    return t;
}

std::expected<void, base::Error> Builder::models(assemble::PhysicsModels& m) const {
    const Section* s = single("models");
    if (!s) return {};
    const std::map<std::string_view, bool assemble::PhysicsModels::*> switches{
        {"doping_mobility", &assemble::PhysicsModels::doping_mobility},
        {"srh", &assemble::PhysicsModels::srh},
        {"auger", &assemble::PhysicsModels::auger},
        {"bgn", &assemble::PhysicsModels::bgn},
        {"field_mobility", &assemble::PhysicsModels::field_mobility},
        {"fermi_dirac", &assemble::PhysicsModels::fermi_dirac},
        {"radiative", &assemble::PhysicsModels::radiative},
        {"incomplete_ionization", &assemble::PhysicsModels::incomplete_ionization},
        {"impact_ionization", &assemble::PhysicsModels::impact_ionization},
        {"btbt_local", &assemble::PhysicsModels::btbt_local},
        {"electrothermal", &assemble::PhysicsModels::electrothermal},
        {"energy_balance", &assemble::PhysicsModels::energy_balance}};
    for (const Entry& e : s->entries) {
        if (const auto it = switches.find(e.key); it != switches.end()) {
            NITCAD_DECK_TRY(v, flag(e));
            m.*(it->second) = v;
        } else if (e.key == "btbt_nonlocal") {
            if (e.value == "off") m.btbt_nonlocal = assemble::NonlocalTunnelling::off;
            else if (e.value == "kane") m.btbt_nonlocal = assemble::NonlocalTunnelling::kane;
            else if (e.value == "direct_wkb") m.btbt_nonlocal = assemble::NonlocalTunnelling::direct_wkb;
            else return std::unexpected(at(e.line, "btbt_nonlocal is off, kane or direct_wkb"));
        } else if (e.key == "hole_mobility") {
            if (e.value == "caughey_thomas") m.hole_mobility = assemble::HoleMobility::caughey_thomas;
            else if (e.value == "masetti") m.hole_mobility = assemble::HoleMobility::masetti;
            else if (e.value == "automatic") m.hole_mobility = assemble::HoleMobility::automatic;
            else return std::unexpected(at(e.line, "hole_mobility is automatic, caughey_thomas or masetti"));
        } else if (e.key == "electron_mobility") {
            if (e.value == "caughey_thomas") m.electron_mobility = assemble::ElectronMobility::caughey_thomas;
            else if (e.value == "masetti") m.electron_mobility = assemble::ElectronMobility::masetti;
            else return std::unexpected(at(e.line, "electron_mobility is caughey_thomas or masetti"));
        } else if (e.key == "surface_mobility") {
            if (e.value == "off") m.surface_mobility = assemble::SurfaceMobility::off;
            else if (e.value == "lombardi") m.surface_mobility = assemble::SurfaceMobility::lombardi;
            else return std::unexpected(at(e.line, "surface_mobility is off or lombardi"));
        } else if (e.key == "quantum_confinement") {
            if (e.value == "off") m.quantum_confinement = assemble::QuantumConfinement::off;
            else if (e.value == "mlda") m.quantum_confinement = assemble::QuantumConfinement::mlda;
            else return std::unexpected(at(e.line, "quantum_confinement is off or mlda"));
        } else if (e.key == "mlda_mass_n") {
            NITCAD_DECK_TRY(v, num(e));
            m.mlda_mass_n = v;
        } else if (e.key == "mlda_mass_p") {
            NITCAD_DECK_TRY(v, num(e));
            m.mlda_mass_p = v;
        } else if (e.key == "ionization_levels") {
            if (e.value == "fixed") m.ionization_levels = assemble::IonizationLevels::fixed;
            else if (e.value == "altermatt") m.ionization_levels = assemble::IonizationLevels::altermatt;
            else return std::unexpected(at(e.line, "ionization_levels is fixed or altermatt"));
        } else if (e.key == "impact_current_resolution") {
            NITCAD_DECK_TRY(v, num(e));
            m.impact_current_resolution = v;
        } else {
            return std::unexpected(unknown(e, *s));
        }
    }
    return {};
}

std::expected<void, base::Error> Builder::solver(Plan& plan) const {
    const Section* s = single("solver");
    if (!s) return {};
    solve::BiasOptions& o = plan.options;
    for (const Entry& e : s->entries) {
        if (e.key == "equations") {
            if (e.value == "drift_diffusion") o.equations = solve::Equations::drift_diffusion;
            else if (e.value == "equilibrium_poisson") o.equations = solve::Equations::equilibrium_poisson;
            else return std::unexpected(at(e.line, "equations is drift_diffusion or equilibrium_poisson"));
        } else if (e.key == "max_iterations") {
            NITCAD_DECK_TRY(v, count(e, 1, 100000));
            o.newton.max_iterations = v;
        } else if (e.key == "tol_update") {
            NITCAD_DECK_TRY(v, num(e));
            o.newton.tol_update = v;
        } else if (e.key == "max_update") {
            NITCAD_DECK_TRY(v, num(e));
            o.newton.max_update = v;
        } else if (e.key == "line_search") {
            if (e.value == "off") o.newton.line_search.enabled = false;
            else if (e.value == "residual" || e.value == "natural") {
                o.newton.line_search.enabled = true;
                o.newton.line_search.test = e.value == "residual" ? solve::LineSearchTest::residual
                                                                  : solve::LineSearchTest::natural;
            } else {
                return std::unexpected(at(e.line, "line_search is off, residual or natural"));
            }
        } else if (e.key == "line_search.max_backtracks") {
            NITCAD_DECK_TRY(v, count(e, 0, 60));
            o.newton.line_search.max_backtracks = v;
        } else if (e.key == "line_search.sufficient_decrease") {
            NITCAD_DECK_TRY(v, num(e));
            o.newton.line_search.sufficient_decrease = v;
        } else if (e.key == "max_cuts") {
            NITCAD_DECK_TRY(v, count(e, 0, 52));
            o.step_control.max_cuts = v;
        } else if (e.key == "growth") {
            NITCAD_DECK_TRY(v, num(e));
            o.step_control.growth = v;
        } else if (e.key == "cut_max_iterations") {
            NITCAD_DECK_TRY(v, count(e, 1, 100000));
            o.step_control.max_iterations = v;
        } else if (e.key == "backend") {
            if (e.value == "eigen") o.linear.backend = linalg::SolverBackend::eigen_sparse_lu;
            else if (e.value == "pardiso") o.linear.backend = linalg::SolverBackend::mkl_pardiso;
            else if (e.value == "gmres_ilut") o.linear.backend = linalg::SolverBackend::gmres_ilut;
            else return std::unexpected(at(e.line, "backend is eigen, pardiso or gmres_ilut"));
        } else if (e.key == "restart") {
            NITCAD_DECK_TRY(v, count(e, 1, 1000));
            o.linear.restart = v;
        } else if (e.key == "max_linear_iterations") {
            NITCAD_DECK_TRY(v, count(e, 1, 1'000'000'000));
            o.linear.max_iterations = v;
        } else if (e.key == "ilut_drop") {
            NITCAD_DECK_TRY(v, num(e));
            if (!(v >= 0.0 && v < 1.0)) return std::unexpected(at(e.line, "ilut_drop is in [0, 1)"));
            o.linear.ilut_drop_tolerance = v;
        } else if (e.key == "ilut_fill") {
            NITCAD_DECK_TRY(v, count(e, 1, 1000));
            o.linear.ilut_fill_factor = v;
        } else if (e.key == "linear_fallback") {
            if (e.value == "sparse_lu") o.linear.fallback = linalg::LinearFallback::sparse_lu;
            else if (e.value == "none") o.linear.fallback = linalg::LinearFallback::none;
            else return std::unexpected(at(e.line, "linear_fallback is sparse_lu or none"));
        } else if (e.key == "mkl") {
            plan.mkl_runtime = file(e);
        } else if (e.key == "threads") {
            NITCAD_DECK_TRY(v, count(e, 1, 1024));
            o.linear.threads = v;
        } else if (e.key == "equilibrate") {
            NITCAD_DECK_TRY(v, flag(e));
            o.linear.equilibrate = v;
        } else if (e.key == "max_backward_error") {
            NITCAD_DECK_TRY(v, num(e));
            o.linear.max_backward_error = v;
        } else if (e.key == "max_refinement_steps") {
            NITCAD_DECK_TRY(v, count(e, 0, 100));
            o.linear.max_refinement_steps = v;
        } else if (e.key == "min_pivot_ratio") {
            NITCAD_DECK_TRY(v, num(e));
            o.linear.min_pivot_ratio = v;
        } else if (e.key == "Ns") {
            NITCAD_DECK_TRY(v, num(e));
            o.Ns_override = v;
        } else if (e.key == "weighted_current") {
            NITCAD_DECK_TRY(v, flag(e));
            o.weighted_current = v;
        } else if (e.key == "minority_floor") {
            NITCAD_DECK_TRY(v, flag(e));
            o.minority_floor = v;
        } else if (e.key == "rates") {
            NITCAD_DECK_TRY(v, flag(e));
            o.store_rates = v;
        } else {
            return std::unexpected(unknown(e, *s));
        }
    }
    if (o.linear.backend == linalg::SolverBackend::mkl_pardiso && plan.mkl_runtime.empty()) {
        return std::unexpected(at(s->line, "backend = pardiso needs 'mkl', the path of mkl_rt"));
    }
    return {};
}

std::expected<void, base::Error> Builder::analysis(Plan& plan) const {
    if (const Section* s = single("equilibrium")) {
        if (!s->entries.empty()) return std::unexpected(unknown(s->entries.front(), *s));
        plan.analysis = Analysis::equilibrium;
        return {};
    }
    if (const Section* s = single("transient")) return transient(*s, plan);
    if (const Section* s = single("small_signal")) return small_signal(*s, plan);
    if (const Section* s = single("trace")) return trace(*s, plan);
    const Section& s = *single("sweep");
    plan.analysis = Analysis::sweep;
    return operating_points(s, plan, [&](const Entry& e) -> std::expected<bool, base::Error> {
        if (e.key != "restart") return false;
        plan.restart = file(e);
        return true;
    });
}

// A refusal of `restart` outside [sweep] (DK9, Unit 43 AN6).
base::Error Builder::no_restart(const Entry& e, const Section& s) const {
    return at(e.line, "restart is only for [sweep] (DK9), not [" + s.kind + "]");
}

std::expected<std::size_t, base::Error> Builder::contact_index(const Entry& e, std::string_view name,
                                                               const Plan& plan) const {
    const auto contacts = plan.device.contacts();
    const auto it = std::find_if(contacts.begin(), contacts.end(), [&](const auto& c) { return c.name == name; });
    if (it == contacts.end()) return std::unexpected(at(e.line, "no contact named '" + std::string(name) + "'"));
    return static_cast<std::size_t>(it - contacts.begin());
}

std::expected<std::size_t, base::Error> Builder::thermal_index(const Entry& e, std::string_view name,
                                                               const Plan& plan) const {
    const auto thermal = plan.device.thermal_contacts();
    const auto it = std::find_if(thermal.begin(), thermal.end(), [&](const auto& t) { return t.name == name; });
    if (it == thermal.end()) return std::unexpected(at(e.line, "no thermal contact named '" + std::string(name) + "'"));
    return static_cast<std::size_t>(it - thermal.begin());
}

// `thermal <name> = <K>` entries of a run that holds one temperature list (a transient, a trace):
// the thermal contacts' temperatures, each contact's own when not named; set only when one is
// named. Returns whether `e` was such an entry.
std::expected<bool, base::Error> Builder::held_temperature(const Entry& e, Plan& plan,
                                                           std::vector<bool>& given) const {
    const auto w = tokens(e.key);
    if (w.size() != 2 || w[0] != "thermal") return false;
    NITCAD_DECK_TRY(index, thermal_index(e, w[1], plan));
    if (given.empty()) given.assign(plan.device.thermal_contacts().size(), false);
    if (given[index]) return std::unexpected(at(e.line, "a second temperature for '" + std::string(w[1]) + "'"));
    given[index] = true;
    NITCAD_DECK_TRY(v, num(e));
    if (plan.options.thermal_bias_K.empty()) {
        std::vector<double> t;
        for (const device::ThermalContact& c : plan.device.thermal_contacts()) t.push_back(c.temperature_K);
        plan.options.thermal_bias_K.push_back(std::move(t));
    }
    plan.options.thermal_bias_K.front()[index] = v;
    return true;
}

std::expected<void, base::Error> Builder::transient(const Section& s, Plan& plan) const {
    plan.analysis = Analysis::transient;
    solve::TransientOptions& o = plan.transient;
    const auto contacts = plan.device.contacts();
    std::vector<std::optional<solve::Waveform>> waveforms(contacts.size());
    std::vector<bool> given_t;
    bool t_end = false;
    for (const Entry& e : s.entries) {
        NITCAD_DECK_TRY(thermal, held_temperature(e, plan, given_t));
        if (thermal) continue;
        const auto w = tokens(e.key);
        if (w.size() == 2 && w[0] == "contact") {
            NITCAD_DECK_TRY(index, contact_index(e, w[1], plan));
            if (waveforms[index]) return std::unexpected(at(e.line, "a second waveform for '" + std::string(w[1]) + "'"));
            NITCAD_DECK_TRY(f, lift(e.line, waveform(e.value)));
            waveforms[index] = std::move(f);
        } else if (e.key == "t_end") {
            NITCAD_DECK_TRY(v, num(e));
            o.t_end_s = v;
            t_end = true;
        } else if (e.key == "dt_initial" || e.key == "dt_min" || e.key == "dt_max") {
            NITCAD_DECK_TRY(v, num(e));
            (e.key == "dt_initial" ? o.dt_initial_s : e.key == "dt_min" ? o.dt_min_s : o.dt_max_s) = v;
        } else if (e.key == "integrator") {
            if (e.value == "bdf2") o.integrator = solve::Integrator::bdf2;
            else if (e.value == "backward_euler") o.integrator = solve::Integrator::backward_euler;
            else return std::unexpected(at(e.line, "integrator is bdf2 or backward_euler"));
        } else if (e.key == "adaptive" || e.key == "fields_every_step") {
            NITCAD_DECK_TRY(v, lift(e.line, boolean(e.value)));
            (e.key == "adaptive" ? o.adaptive : o.fields_every_step) = v;
        } else if (e.key == "rtol") {
            NITCAD_DECK_TRY(v, num(e));
            o.rtol = v;
        } else if (e.key == "density_ref") {
            NITCAD_DECK_TRY(v, num(e));
            o.density_ref_cm3 = v;
        } else if (e.key == "output_times") {
            NITCAD_DECK_TRY(v, lift(e.line, series(e.value, [](std::string_view t) { return number(t); })));
            o.output_times_s = std::move(v);
        } else if (e.key == "restart") {
            return std::unexpected(no_restart(e, s));
        } else {
            return std::unexpected(unknown(e, s));
        }
    }
    if (!t_end) return std::unexpected(at(s.line, "[transient] needs t_end"));
    for (auto& f : waveforms) plan.waveforms.push_back(f ? std::move(*f) : solve::Waveform::constant(0.0));
    return {};
}

std::expected<void, base::Error> Builder::small_signal(const Section& s, Plan& plan) const {
    plan.analysis = Analysis::small_signal;
    bool frequencies = false;
    NITCAD_DECK_CHECK(operating_points(s, plan, [&](const Entry& e) -> std::expected<bool, base::Error> {
        if (e.key == "restart") return std::unexpected(no_restart(e, s));
        if (e.key != "frequencies" && e.key != "field_frequencies") return false;
        NITCAD_DECK_TRY(v, lift(e.line, series(e.value, [](std::string_view t) { return number(t); })));
        if (e.key == "frequencies") {
            plan.small_signal.frequencies_Hz = std::move(v);
            frequencies = true;
        } else {
            plan.small_signal.field_frequencies_Hz = std::move(v);
        }
        return true;
    }));
    if (!frequencies) return std::unexpected(at(s.line, "[small_signal] needs frequencies"));
    return {};
}

std::expected<void, base::Error> Builder::trace(const Section& s, Plan& plan) const {
    plan.analysis = Analysis::trace;
    solve::TraceOptions& o = plan.trace;
    o.start_V.assign(plan.device.contacts().size(), 0.0);
    std::vector<bool> given(o.start_V.size(), false), given_t;
    bool has_contact = false, has_end = false;
    for (const Entry& e : s.entries) {
        NITCAD_DECK_TRY(thermal, held_temperature(e, plan, given_t));
        if (thermal) continue;
        const auto w = tokens(e.key);
        if (w.size() == 2 && w[0] == "start") {
            NITCAD_DECK_TRY(index, contact_index(e, w[1], plan));
            if (given[index]) return std::unexpected(at(e.line, "a second start for '" + std::string(w[1]) + "'"));
            given[index] = true;
            NITCAD_DECK_TRY(v, num(e));
            o.start_V[index] = v;
        } else if (e.key == "contact") {
            NITCAD_DECK_TRY(index, contact_index(e, trimmed(e.value), plan));
            o.contact = index;
            has_contact = true;
        } else if (e.key == "end") {
            NITCAD_DECK_TRY(v, num(e));
            o.end_V = v;
            has_end = true;
        } else if (e.key == "step" || e.key == "max_step" || e.key == "min_step" ||
                   e.key == "density_floor") {
            NITCAD_DECK_TRY(v, num(e));
            (e.key == "step"       ? o.step_V
             : e.key == "max_step" ? o.max_step_V
             : e.key == "min_step" ? o.min_step_V
                                   : o.density_floor_cm3) = v;
        } else if (e.key == "current_limit") {
            NITCAD_DECK_TRY(v, num(e));
            o.current_limit = v;
        } else if (e.key == "max_points") {
            NITCAD_DECK_TRY(v, count(e, 1, 100'000'000));
            o.max_points = static_cast<std::size_t>(v);
        } else if (e.key == "restart") {
            return std::unexpected(no_restart(e, s));
        } else {
            return std::unexpected(unknown(e, s));
        }
    }
    if (!has_contact) return std::unexpected(at(s.line, "[trace] needs contact, the swept contact"));
    if (!has_end) return std::unexpected(at(s.line, "[trace] needs end, the swept contact's last bias"));
    return {};
}

// The operating points of [sweep] and [small_signal] (DK6): `contact <name>` and `thermal <name>`
// series. `other` handles the section's other keys, returning whether it took the entry.
template <class Other>
std::expected<void, base::Error> Builder::operating_points(const Section& s, Plan& plan, Other&& other) const {
    const auto contacts = plan.device.contacts();
    const auto thermal = plan.device.thermal_contacts();
    std::vector<std::vector<double>> bias(contacts.size(), std::vector<double>{0.0});
    std::vector<std::vector<double>> temperature;
    for (const device::ThermalContact& t : thermal) temperature.push_back({t.temperature_K});
    std::vector<bool> given(contacts.size(), false), given_t(thermal.size(), false);
    bool any_thermal = false;
    std::size_t n = 1;
    int n_line = s.line;
    for (const Entry& e : s.entries) {
        NITCAD_DECK_TRY(taken, other(e));
        if (taken) continue;
        const auto w = tokens(e.key);
        if (w.size() != 2 || (w[0] != "contact" && w[0] != "thermal")) return std::unexpected(unknown(e, s));
        const bool is_thermal = w[0] == "thermal";
        std::size_t index = 0;
        if (is_thermal) {
            const auto it = std::find_if(thermal.begin(), thermal.end(), [&](const auto& t) { return t.name == w[1]; });
            if (it == thermal.end()) return std::unexpected(at(e.line, "no thermal contact named '" + std::string(w[1]) + "'"));
            index = static_cast<std::size_t>(it - thermal.begin());
            if (given_t[index]) return std::unexpected(at(e.line, "a second series for '" + std::string(w[1]) + "'"));
            given_t[index] = true;
            any_thermal = true;
        } else {
            const auto it = std::find_if(contacts.begin(), contacts.end(), [&](const auto& c) { return c.name == w[1]; });
            if (it == contacts.end()) return std::unexpected(at(e.line, "no contact named '" + std::string(w[1]) + "'"));
            index = static_cast<std::size_t>(it - contacts.begin());
            if (given[index]) return std::unexpected(at(e.line, "a second series for '" + std::string(w[1]) + "'"));
            given[index] = true;
        }
        NITCAD_DECK_TRY(v, lift(e.line, series(e.value, [](std::string_view t) { return number(t); })));
        if (v.size() > 1) {
            if (n > 1 && v.size() != n) {
                return std::unexpected(at(e.line, "this series has " + std::to_string(v.size()) +
                                                      " values, the one on line " + std::to_string(n_line) +
                                                      " has " + std::to_string(n)));
            }
            n = v.size();
            n_line = e.line;
        }
        (is_thermal ? temperature : bias)[index] = std::move(v);
    }
    const auto at_point = [](const std::vector<double>& v, std::size_t k) { return v.size() == 1 ? v[0] : v[k]; };
    for (std::size_t k = 0; k < n; ++k) {
        std::vector<double> p;
        for (const auto& v : bias) p.push_back(at_point(v, k));
        plan.points.push_back(std::move(p));
        if (any_thermal) {
            std::vector<double> t;
            for (const auto& v : temperature) t.push_back(at_point(v, k));
            plan.options.thermal_bias_K.push_back(std::move(t));
        }
    }
    return {};
}

std::expected<void, base::Error> Builder::outputs(Plan& plan) const {
    const Section* s = single("output");
    if (!s) return {};
    for (const Entry& e : s->entries) {
        std::filesystem::path* target = e.key == "sweep"          ? &plan.outputs.sweep
                                        : e.key == "device"       ? &plan.outputs.device
                                        : e.key == "mesh"         ? &plan.outputs.mesh
                                        : e.key == "checkpoint"   ? &plan.outputs.checkpoint
                                        : e.key == "table"        ? &plan.outputs.table
                                        : e.key == "transient"    ? &plan.outputs.transient
                                        : e.key == "small_signal" ? &plan.outputs.small_signal
                                                                  : nullptr;
        if (!target) return std::unexpected(unknown(e, *s));
        *target = file(e);
    }
    if (plan.analysis == Analysis::equilibrium && (!plan.outputs.sweep.empty() || !plan.outputs.table.empty())) {
        return std::unexpected(at(s->line, "an equilibrium run writes no sweep or table"));
    }
    // Each analysis writes only its own results (Unit 43, AN5).
    const Outputs& o = plan.outputs;
    const auto refuse = [&](const char* what, const char* analysis) {
        return at(s->line, std::string(analysis) + " run writes no " + what);
    };
    const bool writes_sweep = plan.analysis == Analysis::sweep || plan.analysis == Analysis::trace;
    const char* name = plan.analysis == Analysis::equilibrium ? "an equilibrium"
                       : plan.analysis == Analysis::sweep     ? "a sweep"
                       : plan.analysis == Analysis::transient ? "a transient"
                       : plan.analysis == Analysis::trace     ? "a trace"
                                                              : "a small-signal";
    if (!o.sweep.empty() && !writes_sweep) return std::unexpected(refuse("sweep", name));
    if (!o.transient.empty() && plan.analysis != Analysis::transient) {
        return std::unexpected(refuse("transient", name));
    }
    if (!o.small_signal.empty() && plan.analysis != Analysis::small_signal) {
        return std::unexpected(refuse("small_signal", name));
    }
    if (!o.checkpoint.empty() && plan.analysis == Analysis::transient) {
        return std::unexpected(at(s->line, "a transient run writes no checkpoint: its end state is not "
                                           "a steady state a sweep can restart from"));
    }
    return {};
}

std::expected<Plan, base::Error> Builder::run() {
    NITCAD_DECK_CHECK(check_sections());
    device::DeviceDescription d{.mesh = *mesh::make_tensor_grid(std::vector<double>{0.0, 1.0}),
                                .temperature_K = 300.0};
    if (const Section* s = single("device")) {
        for (const Entry& e : s->entries) {
            if (e.key != "temperature") return std::unexpected(unknown(e, *s));
            NITCAD_DECK_TRY(v, num(e));
            d.temperature_K = v;
        }
    }
    for (const Section& s : d_.sections) {
        if (s.kind != "region") continue;
        NITCAD_DECK_TRY(spec, region(s));
        regions_.push_back(std::move(spec));
    }
    if (regions_.empty()) return std::unexpected(at(1, "the deck declares no [region]"));
    for (const RegionSpec& g : regions_) d.regions.push_back(*g.region);
    NITCAD_DECK_CHECK(build_mesh(d));
    NITCAD_DECK_CHECK(doping(d));
    for (const Section& s : d_.sections) {
        if (s.kind == "contact") {
            NITCAD_DECK_TRY(c, contact(s, d));
            d.contacts.push_back(std::move(c));
        } else if (s.kind == "interface") {
            NITCAD_DECK_TRY(f, interface(s));
            d.interfaces.push_back(std::move(f));
        } else if (s.kind == "thermal_contact") {
            NITCAD_DECK_TRY(t, thermal_contact(s));
            d.thermal_contacts.push_back(std::move(t));
        }
    }
    auto device = device::Device::create(std::move(d));
    if (!device) return std::unexpected(located(1, "device", std::move(device.error())));
    Plan plan{.device = std::move(*device), .generated = std::move(generated_)};
    NITCAD_DECK_CHECK(models(plan.options.models));
    NITCAD_DECK_CHECK(solver(plan));
    NITCAD_DECK_CHECK(analysis(plan));
    NITCAD_DECK_CHECK(outputs(plan));
    return plan;
}

}  // namespace

std::expected<Plan, base::Error> build(const Deck& deck) { return Builder(deck).run(); }

}  // namespace NiTCAD::deck

#undef NITCAD_DECK_TRY
#undef NITCAD_DECK_CHECK
