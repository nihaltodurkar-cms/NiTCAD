#include "NiTCAD/deck/build.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <fstream>
#include <iterator>
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
#include "NiTCAD/io/learn_files.hpp"
#include "NiTCAD/mesh/complex.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/heavy_ion.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/materials.hpp"
#include "NiTCAD/physics/polarization.hpp"
#include "NiTCAD/physics/total_dose.hpp"
#include "NiTCAD/solve/adapt_adjoint.hpp"
#include "NiTCAD/solve/adapt_derived.hpp"
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

// Unit 68 (MAD14): an unknown metal preset, with the nearest name.
std::string unknown_metal(const std::string& value) {
    std::vector<std::string_view> names;
    for (const physics::MetalPreset& m : physics::metal_presets()) names.push_back(m.name);
    const std::string_view nearest = physics::nearest_name(value, names);
    return "unknown metal '" + value + "'" + (nearest.empty() ? "" : " (did you mean '" + std::string(nearest) + "'?)");
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
    std::vector<const Entry*> prisms;    // Unit 26l (SD1): layered generate3d meshes
    std::vector<const Entry*> extrudes;  // Unit 53 (FD13): extrude3d meshes
    const Entry* overlap = nullptr;      // Unit 26m (MT4): overlap = replace | fill, layered meshes only
    bool fill = false;
    std::optional<results::MaterialRecord> record;  // Unit 68 (MAD3): when the run records its material
    std::optional<device::Electrode> metal;          // Unit 68 (MAD6): a metal region's electrode
    // Unit 68 (MAQ7a, slice 1): `composition = linear <x_start> <x_end> <axis> <start> <end>`, the
    // material's fraction linear along the axis between the two coordinates (held at the ends
    // outside them); the region's parameter overrides apply to every layer.
    struct Grading {
        const physics::MaterialEntry* entry = nullptr;
        const Entry* line = nullptr;
        double x_start = 0.0, x_end = 0.0, start = 0.0, end = 0.0;
        int axis = 0;
        std::vector<const Entry*> overrides;
    };
    std::optional<Grading> grading;
    const Entry* relaxation = nullptr;  // Unit 68 slice 2 (MAD9): a nitride layer's strain relaxation
};

// A [material] section built (Unit 68, MAD4): its set and its base's composition.
struct UserMaterial {
    physics::MaterialParameters parameters;
    std::vector<double> composition;
    std::string base;  // the registry base's name
};

class Builder {
public:
    Builder(const Deck& d, const solve::MeshRefinement* refinement) : d_(d), refine_(refinement) {}

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
    // A material parameter by its member path (DK7) set from `e`; true when `e` was one.
    std::expected<bool, base::Error> parameter(physics::MaterialParameters& params, const Entry& e) const;
    // Unit 68 (MAQ7a): each graded region split into one uniform-composition region per node layer.
    std::expected<void, base::Error> grade(device::DeviceDescription& d);
    std::expected<void, base::Error> polarization(device::DeviceDescription& d) const;  // Unit 68 slice 2
    // The registry material `e` names, with its composition arguments (Unit 68, MAD1).
    // `composition`, when given, receives the arguments.
    std::expected<physics::MaterialParameters, base::Error> registry_material(const Entry& e,
                                                                              std::vector<double>* composition = nullptr) const;
    // Unit 68 (MAD4): the [material <name>] sections, inline or from their include file.
    std::expected<void, base::Error> user_materials();
    std::expected<UserMaterial, base::Error> user_material(const Section& s) const;
    std::expected<RegionSpec, base::Error> region(const Section& s) const;
    std::expected<std::vector<double>, base::Error> axis(int line, std::string_view value) const;
    std::expected<std::vector<double>, base::Error> axis(const Entry& e) const { return axis(e.line, e.value); }
    std::expected<void, base::Error> build_mesh(device::DeviceDescription& d);
    std::expected<void, base::Error> tensor_mesh(const Section& s, device::DeviceDescription& d);
    std::expected<void, base::Error> generated_2d(const Section& s);
    std::expected<void, base::Error> generated_3d(const Section& s);
    std::expected<void, base::Error> layered_3d(const Section& s);  // Unit 26l (SD1)
    std::expected<void, base::Error> extruded_3d(const Section& s);
    // The meshing keys generate3d and extrude3d share; true when `e` was one of them.
    std::expected<bool, base::Error> option_3d(const Entry& e, mesh::MeshOptions3D& o) const;
    // The swept mesh's keys (Unit 54, UD9, UD10); true when `e` was one of them.
    std::expected<bool, base::Error> option_swept(const Entry& e, mesh::SweptOptions& o) const;
    std::expected<mesh::BoundaryLayer, base::Error> layer(const Entry& e) const;
    std::expected<void, base::Error> file_mesh(const Section& s, device::DeviceDescription& d);
    std::expected<void, base::Error> doping(device::DeviceDescription& d) const;
    std::expected<device::Contact, base::Error> contact(const Section& s,
                                                        const device::DeviceDescription& d) const;
    std::expected<device::Interface, base::Error> interface(const Section& s) const;
    // Unit 69 (TDD2, TDD3): a gate's or an interface's dose keys, gathered in its section's loop.
    struct DoseKeys {
        const Entry* dose = nullptr;
        const Entry* field = nullptr;
        const Entry* trapping = nullptr;
        const Entry* thickness = nullptr;   // interface only
        const Entry* trap_yield = nullptr;  // interface only
    };
    struct DoseCharges {
        double trapped_cm2 = 0.0;  // N_ot
        double traps_cm2 = 0.0;    // Delta N_it
    };
    static bool dose_key(const Entry& e, DoseKeys& k);
    // The charges the keys give (0 and 0 without them). A gate passes its oxide thickness; an
    // interface passes nullopt and takes `dose_oxide_thickness`.
    std::expected<DoseCharges, base::Error> dose_charges(const DoseKeys& k, std::optional<double> gate_thickness_cm) const;
    std::expected<device::ThermalContact, base::Error> thermal_contact(const Section& s) const;
    std::expected<void, base::Error> models(assemble::PhysicsModels& m) const;
    std::expected<void, base::Error> solver(Plan& plan) const;
    std::expected<void, base::Error> analysis(Plan& plan) const;
    std::expected<void, base::Error> transient(const Section& s, Plan& plan) const;
    std::expected<void, base::Error> heavy_ions(Plan& plan) const;  // Unit 65
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
    std::expected<void, base::Error> starts(const Entry* initial, const Entry* restart, Plan& plan) const;
    std::expected<void, base::Error> outputs(Plan& plan) const;
    std::expected<void, base::Error> adapt(Plan& plan) const;
    std::expected<Box, base::Error> box(const Entry& e, int dim) const;

    const Deck& d_;
    std::map<std::string, UserMaterial, std::less<>> user_materials_;
    std::vector<RegionSpec> regions_;
    std::vector<results::MaterialRecord> graded_records_;  // Unit 68: the layers a grading adds
    std::vector<std::size_t> graded_parent_;  // the regions_ index of each graded_records_ layer
    MeshType type_ = MeshType::tensor;
    std::string type_name_;  // the [mesh] type as written (extrude3d for a file mesh from a complex)
    bool layered_mesh_ = false;  // a layered generate2d or generate3d mesh (the only ones taking a region's overlap)
    std::optional<mesh::GeneratedMesh> generated_;
    // Unit 17: the loop's refinement (null: the deck's mesh) and what [adapt] needs of the mesh.
    const solve::MeshRefinement* refine_ = nullptr;
    AdaptPlan mesh_info_;
    bool swept_ = false;
};

std::expected<void, base::Error> Builder::check_sections() const {
    static const std::set<std::string_view> once{"device",       "mesh",      "doping",
                                                 "models",       "solver",    "equilibrium",
                                                 "sweep",        "transient", "small_signal",
                                                 "trace",        "output",    "adapt"};
    static const std::set<std::string_view> many{"region", "contact", "interface", "heavy_ion",
                                                 "thermal_contact", "material"};
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

std::expected<bool, base::Error> Builder::parameter(physics::MaterialParameters& params, const Entry& e) const {
    if (auto* sp = std::get_if<physics::SemiconductorParameters>(&params)) {
        if (e.key == "band_to_band.direct_gap") {
            NITCAD_DECK_TRY(v, flag(e));
            sp->band_to_band.direct_gap = v;
            return true;
        }
        if (double* target = physics::semiconductor_parameter(*sp, e.key)) {
            NITCAD_DECK_TRY(v, num(e));
            *target = v;
            return true;
        }
        return false;
    }
    if (double* target = physics::insulator_parameter(std::get<physics::InsulatorParameters>(params), e.key)) {
        NITCAD_DECK_TRY(v, num(e));
        *target = v;
        return true;
    }
    return false;
}

std::expected<physics::MaterialParameters, base::Error> Builder::registry_material(const Entry& e,
                                                                                   std::vector<double>* composition) const {
    const auto words = tokens(e.value);
    const std::string_view kind = words.front();
    const physics::MaterialEntry* entry = physics::find_material(kind, words.size() - 1);
    if (!entry) {
        // MAC13: the nearest registry or [material] name, when it looks like a misspelling.
        std::vector<std::string_view> names;
        for (const physics::MaterialEntry& m : physics::material_registry()) names.push_back(m.name);
        for (const auto& [name, user] : user_materials_) names.push_back(name);
        const std::string_view nearest = physics::nearest_name(kind, names);
        return std::unexpected(at(e.line, "unknown material '" + std::string(kind) + "'" +
                                              (nearest.empty() ? "" : " (did you mean '" + std::string(nearest) + "'?)")));
    }
    const std::size_t count = words.size() - 1;
    if (count < entry->arguments_min || count > entry->arguments_max) {
        if (entry->arguments_max == 0) return std::unexpected(at(e.line, "'" + std::string(kind) + "' takes no arguments"));
        return std::unexpected(at(e.line, std::string(entry->usage)));
    }
    std::vector<double> arguments;
    for (std::size_t i = 1; i < words.size(); ++i) {
        NITCAD_DECK_TRY(v, lift(e.line, number(words[i])));
        arguments.push_back(v);
    }
    auto made = entry->make(arguments);
    if (!made) return std::unexpected(located(e.line, "material", std::move(made.error())));
    if (composition) *composition = std::move(arguments);
    return std::move(*made);
}

// A [material] section's own entries: base = <registry material> and member-path overrides,
// validated as a region's material is (the error names the parameter).
std::expected<UserMaterial, base::Error> Builder::user_material(const Section& s) const {
    const Entry* base = nullptr;
    for (const Entry& e : s.entries) {
        if (e.key != "base") continue;
        if (base) return std::unexpected(at(e.line, "a second 'base'"));
        base = &e;
    }
    if (!base) return std::unexpected(at(s.line, "[material " + s.name + "] has no 'base'"));
    std::vector<double> composition;
    NITCAD_DECK_TRY(params, registry_material(*base, &composition));
    if (const auto* i = std::get_if<physics::InsulatorParameters>(&params); i && i->metal) {
        return std::unexpected(at(base->line, "'metal' is not a [material] base: a region names it directly"));
    }
    for (const Entry& e : s.entries) {
        if (&e == base) continue;
        if (e.key == "include") return std::unexpected(at(e.line, "an included material file cannot include another"));
        NITCAD_DECK_TRY(taken, parameter(params, e));
        if (!taken) {
            return std::unexpected(at(e.line, "unknown key '" + e.key + "' in [material " + s.name +
                                                  "] (not a parameter of its base)"));
        }
    }
    const std::string where = "[material " + s.name + "]";
    if (const auto* sp = std::get_if<physics::SemiconductorParameters>(&params)) {
        if (auto m = physics::Semiconductor::create(*sp); !m) {
            return std::unexpected(located(base->line, where, std::move(m.error())));
        }
    } else if (auto m = physics::Insulator::create(std::get<physics::InsulatorParameters>(params)); !m) {
        return std::unexpected(located(base->line, where, std::move(m.error())));
    }
    return UserMaterial{std::move(params), std::move(composition), std::string(tokens(base->value).front())};
}

std::expected<void, base::Error> Builder::user_materials() {
    for (const Section& s : d_.sections) {
        if (s.kind != "material") continue;
        if (physics::find_material(s.name)) {
            return std::unexpected(at(s.line, "[material " + s.name + "]: '" + s.name +
                                                  "' is a registry material; a user material needs a name of its own"));
        }
        if (user_materials_.contains(s.name)) return std::unexpected(at(s.line, "a second [material " + s.name + "]"));
        const Entry* include = nullptr;
        for (const Entry& e : s.entries) {
            if (e.key == "include") include = &e;
        }
        if (!include) {
            NITCAD_DECK_TRY(params, user_material(s));
            user_materials_.emplace(s.name, std::move(params));
            continue;
        }
        if (s.entries.size() != 1) {
            return std::unexpected(at(include->line, "a [material] with 'include' takes no other key"));
        }
        // The file holds the section's key lines; parsed under the section's header, its lines
        // renumbered to the file's own, and built by a builder whose messages name the file.
        const std::filesystem::path path = file(*include);
        std::ifstream in(path, std::ios::binary);
        if (!in) return std::unexpected(at(include->line, "cannot open the material file '" + path.string() + "'"));
        const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        std::string_view body = text;
        if (body.starts_with("\xEF\xBB\xBF")) body.remove_prefix(3);
        auto parsed = parse("[material " + s.name + "]\n" + std::string(body), path.filename().string(),
                            path.parent_path());
        if (!parsed) return std::unexpected(located(include->line, "include", std::move(parsed.error())));
        if (parsed->sections.size() != 1) {
            return std::unexpected(at(include->line, "include: a material file holds key lines only, no section"));
        }
        for (Entry& e : parsed->sections.front().entries) --e.line;
        NITCAD_DECK_TRY(params, Builder(*parsed, nullptr).user_material(parsed->sections.front()));
        user_materials_.emplace(s.name, std::move(params));
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
    physics::MaterialParameters params;
    const auto words = tokens(material->value);
    const Entry* composition = nullptr;
    for (const Entry& e : s.entries) {
        if (e.key != "composition") continue;
        if (composition) return std::unexpected(at(e.line, "a second 'composition'"));
        composition = &e;
    }
    if (composition) {  // Unit 68 (MAQ7a): a graded region
        const physics::MaterialEntry* entry = physics::find_material(words.front(), 1);
        if (!entry || entry->arguments_min != 1 || words.size() != 1) {
            return std::unexpected(at(composition->line, "a graded region names a material that takes a fraction, written "
                                                      "without it (for example 'material = sige_strained')"));
        }
        const auto t = tokens(composition->value);
        if (t.size() != 6 || t[0] != "linear") {
            return std::unexpected(at(composition->line, "composition = linear <x_start> <x_end> <axis> <start> <end>"));
        }
        RegionSpec::Grading g{.entry = entry, .line = composition};
        NITCAD_DECK_TRY(x0, lift(composition->line, number(t[1])));
        NITCAD_DECK_TRY(x1, lift(composition->line, number(t[2])));
        if (t[3] != "x" && t[3] != "y" && t[3] != "z") return std::unexpected(at(composition->line, "the grading axis is x, y or z"));
        NITCAD_DECK_TRY(a, lift(composition->line, length(t[4])));
        NITCAD_DECK_TRY(b, lift(composition->line, length(t[5])));
        if (!(b > a)) return std::unexpected(at(composition->line, "the grading's end must lie beyond its start"));
        g.x_start = x0;
        g.x_end = x1;
        g.axis = t[3][0] - 'x';
        g.start = a;
        g.end = b;
        for (const double x : {x0, x1}) {  // both ends must be valid compositions
            const std::array<double, 1> arg{x};
            if (auto made = entry->make(arg); !made) return std::unexpected(located(composition->line, "composition", std::move(made.error())));
        }
        const std::array<double, 1> arg{x0};
        params = *entry->make(arg);
        spec.record = results::MaterialRecord{s.name, std::string(entry->name), {x0}, 0};
        spec.grading = std::move(g);
    } else if (const auto user = user_materials_.find(words.front()); user != user_materials_.end()) {
        if (words.size() != 1) {
            return std::unexpected(at(material->line, "'" + std::string(words.front()) + "' takes no arguments"));
        }
        params = user->second.parameters;
        spec.record = results::MaterialRecord{s.name, user->first, user->second.composition, 0};
    } else {
        std::vector<double> arguments;
        NITCAD_DECK_TRY(made, registry_material(*material, &arguments));
        params = std::move(made);
        if (physics::find_material(words.front(), arguments.size())->recorded) {
            spec.record = results::MaterialRecord{s.name, std::string(words.front()), std::move(arguments), 0};
        }
    }
    const bool is_metal = std::holds_alternative<physics::InsulatorParameters>(params) &&
                          std::get<physics::InsulatorParameters>(params).metal;
    for (const Entry& e : s.entries) {
        if (&e == material || &e == composition) continue;
        if (e.key == "box") { spec.boxes.push_back(&e); continue; }
        if (e.key == "polygon") { spec.polygons.push_back(&e); continue; }
        if (e.key == "seed") { spec.seeds.push_back(&e); continue; }
        if (e.key == "extrude") { spec.extrudes.push_back(&e); continue; }
        if (e.key == "prism") { spec.prisms.push_back(&e); continue; }
        if (e.key == "relaxation") {  // Unit 68 slice 2 (MAD9), checked with polarization
            if (spec.relaxation) return std::unexpected(at(e.line, "a second 'relaxation'"));
            spec.relaxation = &e;
            continue;
        }
        if (e.key == "overlap") {
            if (spec.overlap) return std::unexpected(at(e.line, "a second 'overlap'"));
            if (e.value != "replace" && e.value != "fill") return std::unexpected(at(e.line, "overlap is 'replace' or 'fill'"));
            spec.overlap = &e;
            spec.fill = e.value == "fill";
            continue;
        }
        if (is_metal) {  // Unit 68 (MAD6): the metal's work function, and nothing else
            if (e.key != "work_function" && e.key != "metal") {
                return std::unexpected(at(e.line, "a metal region takes 'work_function' or 'metal', not '" + e.key + "'"));
            }
            if (spec.metal) return std::unexpected(at(e.line, "a metal region takes one of 'work_function' and 'metal'"));
            if (e.key == "work_function") {
                NITCAD_DECK_TRY(v, num(e));
                spec.metal = device::Electrode{.kind = device::GateElectrode::metal, .work_function_eV = v};
                continue;
            }
            const physics::MetalPreset* p = physics::find_metal_preset(e.value);
            if (!p) return std::unexpected(at(e.line, unknown_metal(e.value)));
            if (p->kind == physics::MetalPresetKind::silicide) {
                return std::unexpected(at(e.line, "'" + e.value + "' is a silicide barrier, for a Schottky contact, not a metal region"));
            }
            spec.metal = device::Electrode{
                .kind = p->kind == physics::MetalPresetKind::n_poly   ? device::GateElectrode::n_poly
                        : p->kind == physics::MetalPresetKind::p_poly ? device::GateElectrode::p_poly
                                                                      : device::GateElectrode::metal,
                .work_function_eV = p->work_function_eV};
            continue;
        }
        NITCAD_DECK_TRY(taken, parameter(params, e));
        if (taken) {
            if (spec.grading) spec.grading->overrides.push_back(&e);
            continue;
        }
        if (e.key == "metal" && !physics::find_metal_preset(e.value)) return std::unexpected(at(e.line, unknown_metal(e.value)));
        if (e.key == "work_function" || e.key == "metal") {
            return std::unexpected(at(e.line, "'" + e.key + "' is for a metal region (material = metal)"));
        }
        return std::unexpected(at(e.line, "unknown key '" + e.key + "' in [region " + s.name +
                                              "] (not a parameter of its material)"));
    }
    if (is_metal && !spec.metal) {
        return std::unexpected(at(material->line, "a metal region needs 'work_function' or 'metal = <preset>'"));
    }
    if (spec.record) spec.record->parameter_hash = physics::parameter_hash(params);
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
std::expected<std::vector<double>, base::Error> Builder::axis(int line, std::string_view v) const {
    std::vector<std::string_view> pieces;
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
            if (!piece.ends_with(')')) return std::unexpected(at(line, "linspace( without ')'"));
            const auto args = items(piece.substr(9, piece.size() - 10));
            if (args.size() != 3) return std::unexpected(at(line, "linspace takes (first, last, count)"));
            NITCAD_DECK_TRY(a, lift(line, length(args[0])));
            NITCAD_DECK_TRY(b, lift(line, length(args[1])));
            NITCAD_DECK_TRY(n, lift(line, integer(args[2])));
            if (n < 2 || n > 10'000'000) return std::unexpected(at(line, "linspace needs a count of at least 2"));
            for (long long k = 0; k < n; ++k) {
                append(k == n - 1 ? b : a + (b - a) * static_cast<double>(k) / static_cast<double>(n - 1));
            }
            continue;
        }
        for (const std::string_view item : items(piece)) {
            NITCAD_DECK_TRY(x, lift(line, length(item)));
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
    // Unit 17: an axis the loop refined is taken whole (its straddles already in it); a straddle's
    // two lines and the span between them are protected from refinement either way.
    const auto refined = [&](std::size_t a) { return refine_ != nullptr && !refine_->axes[a].empty(); };
    for (const Entry* e : straddles) {
        const auto which = e->key.substr(9);
        const auto it = std::find(names.begin(), names.begin() + dim, which);
        if (it == names.begin() + dim) return std::unexpected(at(e->line, "straddle on an axis the mesh does not have"));
        NITCAD_DECK_TRY(ps, lengths(*e, e->value, 2));
        const auto k = static_cast<std::size_t>(it - names.begin());
        mesh_info_.axis_protected[k].push_back({ps[0] - 0.5 * ps[1], ps[0] + 0.5 * ps[1]});
        if (refined(k)) continue;
        auto& a = *axes[k];
        auto moved = mesh::straddle_interface(std::move(a), ps[0], ps[1]);
        if (!moved) return std::unexpected(located(e->line, "straddle", std::move(moved.error())));
        a = std::move(*moved);
    }
    for (std::size_t k = 0; k < static_cast<std::size_t>(dim); ++k) {
        if (refined(k)) axes[k] = refine_->axes[k];
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
    bool layered = false;            // Unit 26k (GD1): the regions' shapes layered in deck order
    const Entry* domain = nullptr;   // the layered domain, x0 x1 y0 y1 (default: the shapes' bounds)
    std::vector<const Entry*> patch_entries;  // in the mesher's patch order (its errors name a patch by index)
    for (const Entry& e : s.entries) {
        if (e.key == "type") continue;
        if (e.key == "layered") {
            NITCAD_DECK_TRY(v, flag(e));
            layered = v;
        } else if (e.key == "domain") {
            domain = &e;
        } else if (e.key.starts_with("patch ")) {
            NITCAD_DECK_TRY(ab, points(e, e.value, 2));
            if (ab.size() != 2) return std::unexpected(at(e.line, "a patch is two points"));
            g.patches.push_back({e.key.substr(6), ab[0], ab[1]});
            patch_entries.push_back(&e);
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
    if (domain != nullptr && !layered) return std::unexpected(at(domain->line, "domain is for a layered generate2d mesh (layered = on)"));
    if (layered) {
        // GD1: one region without shapes is the background; every box and polygon, in deck order, replaces what it
        // covers (or, in a region with overlap = fill, takes only what no earlier shape covers: Unit 26m MT4); the
        // kernel resolves the layers into the tiling (GD2).
        layered_mesh_ = true;
        mesh::LayeredGeometry2D l;
        std::vector<std::pair<const Entry*, std::int32_t>> shapes;
        for (std::size_t r = 0; r < regions_.size(); ++r) {
            const RegionSpec& spec = regions_[r];
            if (!spec.seeds.empty()) return std::unexpected(at(spec.seeds.front()->line, "a generate2d mesh takes no seed"));
            if (spec.boxes.empty() && spec.polygons.empty()) {
                if (l.background >= 0) {
                    return std::unexpected(at(spec.section->line, "two regions without a shape: only one may be the background"));
                }
                l.background = static_cast<std::int32_t>(r);
            }
            for (const Entry* e : spec.boxes) shapes.emplace_back(e, static_cast<std::int32_t>(r));
            for (const Entry* e : spec.polygons) shapes.emplace_back(e, static_cast<std::int32_t>(r));
        }
        std::ranges::sort(shapes, {}, [](const auto& p) { return p.first->line; });
        double lo[2] = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
        double hi[2] = {-lo[0], -lo[1]};
        for (const auto& [e, r] : shapes) {
            mesh::LayeredShape shape{{}, r, regions_[static_cast<std::size_t>(r)].fill};
            if (e->key == "box") {
                NITCAD_DECK_TRY(b, box(*e, 2));
                shape.vertices = {{b.bounds[0], b.bounds[2], 0.0}, {b.bounds[1], b.bounds[2], 0.0},
                                  {b.bounds[1], b.bounds[3], 0.0}, {b.bounds[0], b.bounds[3], 0.0}};
            } else {
                NITCAD_DECK_TRY(v, points(*e, e->value, 2));
                shape.vertices = std::move(v);
            }
            for (const auto& p : shape.vertices) {
                for (int k = 0; k < 2; ++k) {
                    lo[k] = std::min(lo[k], p[static_cast<std::size_t>(k)]);
                    hi[k] = std::max(hi[k], p[static_cast<std::size_t>(k)]);
                }
            }
            l.shapes.push_back(std::move(shape));
        }
        if (domain != nullptr) {
            NITCAD_DECK_TRY(b, box(*domain, 2));
            l.lo = {b.bounds[0], b.bounds[2], 0.0};
            l.hi = {b.bounds[1], b.bounds[3], 0.0};
        } else if (l.shapes.empty()) {
            return std::unexpected(at(s.line, "a layered generate2d mesh needs a domain or a shape"));
        } else {
            l.lo = {lo[0], lo[1], 0.0};
            l.hi = {hi[0], hi[1], 0.0};
        }
        l.patches = std::move(g.patches);
        auto tiled = mesh::layered_partition2d(l);
        if (!tiled) {
            const auto& err = tiled.error();
            const int line = err.context && err.context->index && *err.context->index < shapes.size()
                                 ? shapes[*err.context->index].first->line
                                 : s.line;
            return std::unexpected(located(line, "[mesh] layered", std::move(tiled.error())));
        }
        g = std::move(*tiled);
    } else {
        for (std::size_t r = 0; r < regions_.size(); ++r) {
            const RegionSpec& spec = regions_[r];
            if (!spec.boxes.empty() || !spec.seeds.empty()) {
                return std::unexpected(at((spec.boxes.empty() ? spec.seeds : spec.boxes).front()->line,
                                          "a generate2d mesh takes polygon in its regions (box with layered = on)"));
            }
            for (const Entry* e : spec.polygons) {
                NITCAD_DECK_TRY(v, points(*e, e->value, 2));
                g.polygons.push_back({std::move(v), static_cast<std::int32_t>(r)});
            }
        }
    }
    if (refine_ != nullptr) o.size_field = refine_->size_field;
    mesh_info_.kind = solve::AdaptMeshKind::generated2d;
    mesh_info_.interface_resolution = o.interface_resolution;
    mesh_info_.max_size = o.max_size;
    auto m = mesh::generate_mesh(g, o);
    if (!m) {
        // Unit 26k (GD9): an error of one patch (its index the context's) is that patch's line, so the deck pane and
        // the Geometry tab point at it; every other error the section's.
        const auto& err = m.error();
        const bool of_patch = err.message.starts_with("patch") || err.message.starts_with("a patch end within");
        const int line = of_patch && err.context && err.context->index && *err.context->index < patch_entries.size() &&
                                 g.patches.size() == patch_entries.size()
                             ? patch_entries[*err.context->index]->line
                             : s.line;
        return std::unexpected(located(line, "[mesh] generate2d", std::move(m.error())));
    }
    generated_ = std::move(*m);
    return {};
}

std::expected<void, base::Error> Builder::generated_3d(const Section& s) {
    for (const Entry& e : s.entries) {
        if (e.key != "layered") continue;
        NITCAD_DECK_TRY(on, flag(e));
        if (on) return layered_3d(s);
    }
    mesh::Geometry3D g;
    mesh::MeshOptions3D o;
    std::vector<const Entry*> patches;
    bool cuboid = false;
    for (const Entry& e : s.entries) {
        if (e.key == "type" || e.key == "layered") continue;
        if (e.key == "domain") return std::unexpected(at(e.line, "domain is for a layered generate3d mesh (layered = on)"));
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
        if (!spec.boxes.empty() || !spec.polygons.empty() || !spec.prisms.empty()) {
            return std::unexpected(at((!spec.boxes.empty() ? spec.boxes : !spec.polygons.empty() ? spec.polygons : spec.prisms).front()->line,
                                      "a generate3d mesh takes seed in its regions (box and prism with layered = on)"));
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

// Unit 26l (SD1): a layered generate3d mesh. `[mesh]` gives the domain (x0 x1 y0 y1 z0 z1), the patch boxes
// (`patch <name> = x0 y0 z0 x1 y1 z1`, as extrude3d's) and the 3D meshing keys; one region without solids is the
// background, and every region `box` (x0 x1 y0 y1 z0 z1) and `prism` (<axis> <from> <to>, u v, ...) is painted in
// deck order over what it covers; the kernel resolves the layers into the polyhedral complex (SD2).
std::expected<void, base::Error> Builder::layered_3d(const Section& s) {
    layered_mesh_ = true;
    mesh::LayeredGeometry3D l;
    mesh::MeshOptions3D o;
    const Entry* domain = nullptr;
    std::vector<const Entry*> patch_entries;
    for (const Entry& e : s.entries) {
        if (e.key == "type" || e.key == "layered") continue;
        if (e.key == "domain") {
            domain = &e;
        } else if (e.key.starts_with("patch ")) {
            const auto w = tokens(e.key);
            if (w.size() != 2) return std::unexpected(at(e.line, "a layered generate3d patch is 'patch <name> = <x0 y0 z0 x1 y1 z1>'"));
            NITCAD_DECK_TRY(b, lengths(e, e.value, 6));
            if (b[0] > b[3] || b[1] > b[4] || b[2] > b[5]) return std::unexpected(at(e.line, "a box needs lo <= hi on each axis"));
            l.patches.push_back({std::string(w[1]), {b[0], b[1], b[2]}, {b[3], b[4], b[5]}});
            patch_entries.push_back(&e);
        } else if (e.key == "cuboid" || e.key == "vertex" || e.key == "plane" || e.key.starts_with("face ")) {
            return std::unexpected(at(e.line, "'" + e.key + "' is for generate3d without layers: a layered one takes domain"));
        } else {
            NITCAD_DECK_TRY(taken, option_3d(e, o));
            if (!taken) return std::unexpected(unknown(e, s));
        }
    }
    if (domain == nullptr) return std::unexpected(at(s.line, "a layered generate3d mesh needs a domain (x0 x1 y0 y1 z0 z1)"));
    NITCAD_DECK_TRY(d, box(*domain, 3));
    l.lo = {d.bounds[0], d.bounds[2], d.bounds[4]};
    l.hi = {d.bounds[1], d.bounds[3], d.bounds[5]};
    std::vector<std::pair<const Entry*, std::int32_t>> solids;
    for (std::size_t r = 0; r < regions_.size(); ++r) {
        const RegionSpec& spec = regions_[r];
        if (!spec.seeds.empty() || !spec.polygons.empty()) {
            return std::unexpected(at((spec.seeds.empty() ? spec.polygons : spec.seeds).front()->line,
                                      "a layered generate3d mesh takes box and prism in its regions"));
        }
        if (spec.boxes.empty() && spec.prisms.empty()) {
            if (l.background >= 0) return std::unexpected(at(spec.section->line, "two regions without a shape: only one may be the background"));
            l.background = static_cast<std::int32_t>(r);
        }
        for (const Entry* e : spec.boxes) solids.emplace_back(e, static_cast<std::int32_t>(r));
        for (const Entry* e : spec.prisms) solids.emplace_back(e, static_cast<std::int32_t>(r));
    }
    std::ranges::sort(solids, {}, [](const auto& p) { return p.first->line; });
    for (const auto& [e, r] : solids) {
        mesh::LayeredSolid v;
        v.region = r;
        v.fill = regions_[static_cast<std::size_t>(r)].fill;
        if (e->key == "box") {
            NITCAD_DECK_TRY(b, box(*e, 3));
            v.kind = mesh::LayeredSolid::Kind::box;
            v.lo = {b.bounds[0], b.bounds[2], b.bounds[4]};
            v.hi = {b.bounds[1], b.bounds[3], b.bounds[5]};
        } else {
            // prism = <axis> <from> <to>, u v, u v, ...
            const std::string_view text = e->value;
            const std::size_t comma = text.find(',');
            const auto head = tokens(text.substr(0, comma));
            if (comma == std::string_view::npos || head.size() != 3 || (head[0] != "x" && head[0] != "y" && head[0] != "z")) {
                return std::unexpected(at(e->line, "a prism is '<axis x|y|z> <from> <to>, <u v>, <u v>, ...'"));
            }
            v.kind = mesh::LayeredSolid::Kind::prism;
            v.axis = head[0] == "x" ? 0 : head[0] == "y" ? 1 : 2;
            NITCAD_DECK_TRY(from, lift(e->line, length(head[1])));
            NITCAD_DECK_TRY(to, lift(e->line, length(head[2])));
            v.from = from;
            v.to = to;
            NITCAD_DECK_TRY(p, points(*e, text.substr(comma + 1), 2));
            v.polygon = std::move(p);
        }
        l.solids.push_back(std::move(v));
    }
    auto complex = mesh::layered_partition3d(l);
    if (!complex) {
        const auto& err = complex.error();
        const bool of_patch = err.message.starts_with("a patch box");
        const auto& list = of_patch ? patch_entries : std::vector<const Entry*>{};
        int line = s.line;
        if (err.context && err.context->index) {
            const std::size_t k = *err.context->index;
            if (of_patch && k < list.size()) line = list[k]->line;
            else if (!of_patch && k < solids.size()) line = solids[k].first->line;
        }
        return std::unexpected(located(line, "[mesh] layered", std::move(complex.error())));
    }
    auto m = mesh::generate_mesh(*complex, o);
    if (!m) {
        // Unit 26l (SD5): the generator's interface and patch-edge rules name the facet they fail on; a layered
        // deck says where it is (the facet's vertex centroid, cm), so the modeller can mark it.
        auto err = std::move(m.error());
        const std::string_view msg = err.message;
        const bool of_facet = (msg.find("interface") != std::string_view::npos && !msg.starts_with("interface_resolution")) ||
                              msg.starts_with("a patch edge") || msg.starts_with("a patch-edge vertex") || msg.starts_with("a domain-edge vertex");
        if (of_facet && err.context && err.context->index && *err.context->index < complex->facets.size()) {
            const auto& f = complex->facets[*err.context->index];
            std::array<double, 3> c{};
            for (const std::size_t v : f.vertices) {
                for (std::size_t k = 0; k < 3; ++k) c[k] += complex->vertices[v][k] / static_cast<double>(f.vertices.size());
            }
            err.message += std::format(" near ({:.6g}, {:.6g}, {:.6g})", c[0], c[1], c[2]);
        }
        return std::unexpected(located(s.line, "[mesh] generate3d", std::move(err)));
    }
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

std::expected<bool, base::Error> Builder::option_swept(const Entry& e, mesh::SweptOptions& o) const {
    mesh::MeshOptions2D& c = o.section;
    if (e.key == "max_size") {
        NITCAD_DECK_TRY(v, len(e));
        c.max_size = v;
    } else if (e.key == "interface_resolution") {
        NITCAD_DECK_TRY(v, len(e));
        c.interface_resolution = v;
    } else if (e.key == "min_angle") {
        NITCAD_DECK_TRY(v, num(e));
        c.min_angle_deg = v;
    } else if (e.key == "max_nodes") {  // the 3D nodes (UD9)
        NITCAD_DECK_TRY(v, count(e, 1, std::numeric_limits<int>::max()));
        o.max_nodes = static_cast<std::size_t>(v);
    } else if (e.key == "box") {  // in the cross-section: u0 v0 u1 v1 h
        NITCAD_DECK_TRY(b, lengths(e, e.value, 5));
        c.boxes.push_back({{b[0], b[1], 0.0}, {b[2], b[3], 0.0}, b[4]});
    } else if (e.key == "lines") {
        NITCAD_DECK_TRY(v, axis(e));
        o.lines = std::move(v);
    } else if (e.key == "straddle") {  // <station> <width>
        NITCAD_DECK_TRY(b, lengths(e, e.value, 2));
        o.straddles.push_back({b[0], b[1]});
    } else if (e.key == "layer") {
        NITCAD_DECK_TRY(l, layer(e));
        c.layers.push_back(std::move(l));
    } else {
        return false;
    }
    return true;
}

// layer = <u v>, <u v> along=<lengths> lines=<signed lengths> [straddle=<length>] (UD10): the
// base's end points, then the options, each value running to the next option's name; along and
// lines are axes (lengths, linspace pieces joined by '+'), lines signed (positive to the left of
// the base's direction).
std::expected<mesh::BoundaryLayer, base::Error> Builder::layer(const Entry& e) const {
    const std::string_view v = e.value;
    struct Option {
        std::string_view name;
        std::size_t name_at = 0, value_at = 0;
    };
    std::vector<Option> found;
    int depth = 0;
    const auto word = [](char ch) {
        return (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') || (ch >= '0' && ch <= '9') || ch == '_';
    };
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (v[i] == '(') ++depth;
        if (v[i] == ')') --depth;
        if (v[i] != '=' || depth != 0) continue;
        std::size_t k = i;
        while (k > 0 && word(v[k - 1])) --k;
        found.push_back({v.substr(k, i - k), k, i + 1});
    }
    if (found.empty()) {
        return std::unexpected(at(e.line, "a layer is '<u v>, <u v> along=<lengths> lines=<signed lengths> [straddle=<length>]'"));
    }
    mesh::BoundaryLayer l;
    NITCAD_DECK_TRY(ab, points(e, trimmed(v.substr(0, found.front().name_at)), 2));
    if (ab.size() != 2) return std::unexpected(at(e.line, "a layer's base is two points"));
    l.a = ab[0];
    l.b = ab[1];
    bool along = false, lines = false, straddle = false;
    for (std::size_t i = 0; i < found.size(); ++i) {
        const Option& p = found[i];
        const std::size_t end = i + 1 < found.size() ? found[i + 1].name_at : v.size();
        const std::string_view text = trimmed(v.substr(p.value_at, end - p.value_at));
        bool* seen = p.name == "along" ? &along : p.name == "lines" ? &lines : p.name == "straddle" ? &straddle : nullptr;
        if (!seen) {
            return std::unexpected(at(e.line, "unknown layer option '" + std::string(p.name) + "=' (along=, lines=, straddle=)"));
        }
        if (*seen) return std::unexpected(at(e.line, "a second '" + std::string(p.name) + "=' in a layer"));
        *seen = true;
        if (text.empty()) return std::unexpected(at(e.line, "the layer option '" + std::string(p.name) + "=' needs a value"));
        if (p.name == "straddle") {
            NITCAD_DECK_TRY(sv, lift(e.line, length(text)));
            l.straddle = sv;
        } else {
            NITCAD_DECK_TRY(vals, axis(e.line, text));
            (p.name == "along" ? l.along : l.lines) = std::move(vals);
        }
    }
    if (!along || !lines) return std::unexpected(at(e.line, "a layer needs along= and lines="));
    return l;
}

// Unit 53 (DECISIONS.md FD13): a sectioned extrusion. `[mesh]` gives the axis, the stations, the
// cells and the patch boxes; each region assigns its cells per section with `extrude`.
//
// Unit 54 (UD9, UD10): `sweep = on` selects the swept mesh (mesh::generate_swept), whose keys are
// lines, straddle, layer and min_angle, with box in the cross-section's (u, v); radius_edge is the
// Delaunay generator's and refused there, and the swept keys are refused without it.
std::expected<void, base::Error> Builder::extruded_3d(const Section& s) {
    mesh::ExtrudedGeometry x;
    mesh::MeshOptions3D o;
    mesh::SweptOptions sw;
    const Entry* axis_entry = nullptr;
    const Entry* stations_entry = nullptr;
    const Entry* lines_entry = nullptr;
    bool sweep = false;
    for (const Entry& e : s.entries) {
        if (e.key != "sweep") continue;
        NITCAD_DECK_TRY(on, flag(e));
        sweep = on;
    }
    static const std::set<std::string_view> swept_keys{"lines", "straddle", "layer", "min_angle"};
    for (const Entry& e : s.entries) {
        if (e.key == "type" || e.key == "sweep") continue;
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
        } else if (sweep) {
            if (e.key == "radius_edge") {
                return std::unexpected(at(e.line, "radius_edge is for sweep = off: a swept mesh's cross-section is meshed "
                                                  "in 2D (min_angle)"));
            }
            if (e.key == "lines") {
                if (lines_entry) return std::unexpected(at(e.line, "a second 'lines'"));
                lines_entry = &e;
            }
            NITCAD_DECK_TRY(taken, option_swept(e, sw));
            if (!taken) return std::unexpected(unknown(e, s));
        } else {
            if (swept_keys.contains(e.key)) return std::unexpected(at(e.line, "'" + e.key + "' is for a swept mesh (sweep = on)"));
            NITCAD_DECK_TRY(taken, option_3d(e, o));
            if (!taken) return std::unexpected(unknown(e, s));
        }
    }
    if (!axis_entry) return std::unexpected(at(s.line, "an extrude3d mesh needs 'axis'"));
    if (!stations_entry) return std::unexpected(at(s.line, "an extrude3d mesh needs 'stations'"));
    if (x.cells.empty()) return std::unexpected(at(s.line, "an extrude3d mesh needs at least one 'cell'"));
    if (sweep && !lines_entry) return std::unexpected(at(s.line, "a swept extrude3d mesh (sweep = on) needs 'lines'"));
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
            // The nearest station within rounding of the stations' span: 30nm and 0.03um are one
            // position though their doubles may differ in the last place.
            const auto station = [&](double t) -> std::optional<std::size_t> {
                if (x.stations.empty()) return std::nullopt;
                const auto [lo_it, hi_it] = std::ranges::minmax_element(x.stations);
                const double tol = 1e-12 * std::max(*hi_it - *lo_it, std::max(std::abs(*lo_it), std::abs(*hi_it)));
                std::optional<std::size_t> best;
                for (std::size_t i = 0; i < x.stations.size(); ++i) {
                    if (std::abs(x.stations[i] - t) > tol) continue;
                    if (!best || std::abs(x.stations[i] - t) < std::abs(x.stations[*best] - t)) best = i;
                }
                return best;
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
    if (sweep) {
        swept_ = true;
        mesh_info_.kind = solve::AdaptMeshKind::swept;
        mesh_info_.sweep_axis = x.axis;
        for (const mesh::AxisStraddle& st : sw.straddles) {
            mesh_info_.sweep_protected.push_back({st.station - 0.5 * st.width, st.station + 0.5 * st.width});
        }
        mesh_info_.interface_resolution = sw.section.interface_resolution;
        mesh_info_.max_size = sw.section.max_size;
        if (refine_ != nullptr) {
            sw.section.size_field = refine_->size_field;
            if (!refine_->sweep_lines.empty()) sw.lines = refine_->sweep_lines;
        }
        auto m = mesh::generate_swept(x, sw);
        if (!m) return std::unexpected(located(s.line, "[mesh] extrude3d (sweep)", std::move(m.error())));
        generated_ = std::move(*m);
        return {};
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
        // A generation digest belongs to the complex generator, which a deck reaches as extrude3d:
        // the run record names that generator, not "file" (FD14).
        if (generated_->generation_digest != 0) type_name_ = "extrude3d";
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
    type_name_ = type->value;
    if (type_ != MeshType::extrude3d) {
        for (const RegionSpec& g : regions_) {
            if (!g.extrudes.empty()) return std::unexpected(at(g.extrudes.front()->line, "extrude is for an extrude3d mesh"));
        }
    }
    if (type_ != MeshType::generate3d) {
        for (const RegionSpec& g : regions_) {
            if (!g.prisms.empty()) return std::unexpected(at(g.prisms.front()->line, "prism is for a layered generate3d mesh (layered = on)"));
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
    if (!layered_mesh_) {
        for (const RegionSpec& g : regions_) {
            if (g.overlap) return std::unexpected(at(g.overlap->line, "overlap is for a layered mesh (generate2d or generate3d, layered = on)"));
        }
    }
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
            const auto p = parse_profile(e.value);
            if (!p) return std::unexpected(at(e.line, p.error()));
            profiles.push_back(*p);
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


std::expected<void, base::Error> Builder::grade(device::DeviceDescription& d) {
    for (std::size_t r = 0; r < regions_.size(); ++r) {
        RegionSpec& spec = regions_[r];
        if (!spec.grading) continue;
        const RegionSpec::Grading& g = *spec.grading;
        if (g.axis >= d.mesh.dimension()) {
            return std::unexpected(at(g.line->line, std::format("the grading axis lies outside this {}D mesh", d.mesh.dimension())));
        }
        for (const Section& s : d_.sections) {
            if (s.kind == "interface" && s.name.find(spec.section->name) != std::string::npos) {
                for (const std::string_view w : tokens(s.name)) {
                    if (w == spec.section->name) {
                        return std::unexpected(at(s.line, "[interface] cannot name the graded region '" + spec.section->name + "'"));
                    }
                }
            }
        }
        // The node layers: every node of the region by its coordinate along the axis.
        std::map<double, std::vector<std::size_t>> layers;
        const auto points = d.mesh.points();
        for (std::size_t i = 0; i < d.node_region.size(); ++i) {
            if (d.node_region[i] == static_cast<device::RegionId>(r)) layers[points[i][static_cast<std::size_t>(g.axis)]].push_back(i);
        }
        // shortcut: one region per distinct coordinate; an unstructured mesh has one per node, so
        // it is capped; a tensor or extruded mesh has one per grid line.
        constexpr std::size_t max_layers = 512;
        if (layers.size() > max_layers) {
            return std::unexpected(at(g.line->line, std::format("the graded region has {} node layers along its axis (at most {}): "
                                                                "grade on a tensor or extruded mesh",
                                                                layers.size(), max_layers)));
        }
        std::size_t k = 0;
        for (const auto& [c, nodes] : layers) {
            const double t = std::clamp((c - g.start) / (g.end - g.start), 0.0, 1.0);
            const double x = g.x_start + (g.x_end - g.x_start) * t;
            const std::array<double, 1> arg{x};
            auto params = g.entry->make(arg);
            if (!params) return std::unexpected(located(g.line->line, "composition", std::move(params.error())));
            for (const Entry* e : g.overrides) NITCAD_DECK_CHECK(parameter(*params, *e));
            auto made = physics::Semiconductor::create(std::get<physics::SemiconductorParameters>(*params));
            if (!made) return std::unexpected(located(g.line->line, "composition", std::move(made.error())));
            results::MaterialRecord record{spec.section->name, std::string(g.entry->name), {x}, physics::parameter_hash(*params)};
            if (k == 0) {
                d.regions[r].material = std::move(*made);
                spec.record = std::move(record);
            } else {
                const std::string name = spec.section->name + "." + std::to_string(k);
                const auto index = static_cast<device::RegionId>(d.regions.size());
                d.regions.push_back({name, std::move(*made)});
                for (const std::size_t i : nodes) d.node_region[i] = index;
                record.region = name;
                graded_records_.push_back(std::move(record));
                graded_parent_.push_back(r);
            }
            ++k;
        }
    }
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
    const Entry* electrode = nullptr;
    const Entry* metal = nullptr;  // Unit 68 (MAD14): a named metal in place of the work function
    const Entry* region_entry = nullptr;  // Unit 68 (MAD6): the metal region an electrode holds
    const Entry* silicide = nullptr;      // Unit 68 (MAD14, slice 5c): a silicide barrier preset
    bool schottky_keys = false;
    DoseKeys dose;
    for (const Entry& e : s.entries) {
        if (dose_key(e, dose)) gate_keys = true;
        else if (e.key == "kind") kind = &e;
        else if (e.key == "metal") metal = &e;
        else if (e.key == "region") region_entry = &e;
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
            electrode = &e;
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
        return std::unexpected(at(s.line, "oxide_thickness, oxide_permittivity, fixed_charge and the dose keys are for a gate"));
    }
    if (dose.thickness || dose.trap_yield) {
        const Entry& e = dose.thickness ? *dose.thickness : *dose.trap_yield;
        return std::unexpected(at(e.line, "'" + e.key + "' is for an [interface]; a gate's dose uses its oxide_thickness"));
    }
    if (c.kind == device::ContactKind::gate) {
        NITCAD_DECK_TRY(q, dose_charges(dose, c.gate.oxide_thickness_cm));
        if (q.trapped_cm2 != 0.0) c.gate.fixed_charge_cm2 += q.trapped_cm2;  // without a dose the value stays as written (-0 too)
    }
    if ((barrier || schottky_keys) && c.kind != device::ContactKind::schottky) {
        return std::unexpected(at(s.line, "barrier, richardson_n and richardson_p are for a Schottky contact"));
    }
    if (metal) {  // Unit 68 (MAD14): the preset stands for 'electrode' and 'work_function'
        const physics::MetalPreset* p = physics::find_metal_preset(metal->value);
        if (!p) return std::unexpected(at(metal->line, unknown_metal(metal->value)));
        if (work_function || barrier || electrode) {
            return std::unexpected(at(metal->line, "'metal' replaces 'work_function', 'barrier' and 'electrode'; give one"));
        }
        if (c.kind == device::ContactKind::ohmic) {
            return std::unexpected(at(metal->line, "'metal' is for a gate, electrode or Schottky contact"));
        }
        if (p->kind == physics::MetalPresetKind::silicide) {  // slice 5c: a measured barrier (MAD14)
            if (c.kind != device::ContactKind::schottky) {
                return std::unexpected(at(metal->line, "'" + metal->value + "' is a silicide barrier, for a Schottky contact "
                                                                            "on silicon (silicide barriers do not follow Schottky-Mott)"));
            }
            c.schottky.barrier_eV = p->barrier_eV;
            barrier = metal;
            silicide = metal;
        } else if (p->kind != physics::MetalPresetKind::metal && c.kind == device::ContactKind::schottky) {
            return std::unexpected(at(metal->line, "'" + metal->value + "' is a polysilicon electrode, not a Schottky metal"));
        }
        const device::GateElectrode g = p->kind == physics::MetalPresetKind::n_poly   ? device::GateElectrode::n_poly
                                        : p->kind == physics::MetalPresetKind::p_poly ? device::GateElectrode::p_poly
                                                                                      : device::GateElectrode::metal;
        // Exactly as written out: `electrode = <kind>` (gate, electrode) and `work_function = <eV>`.
        if (c.kind != device::ContactKind::schottky) {
            c.gate.electrode = g;
            c.electrode.kind = g;
        }
        if (!silicide) {
            c.gate.work_function_eV = c.electrode.work_function_eV = c.schottky.work_function_eV = p->work_function_eV;
            work_function = metal;  // a Schottky contact takes it as its work function (SD1)
        }
    }
    if (c.kind == device::ContactKind::schottky) {
        // Exactly one of the barrier and the work function (SD1).
        if ((barrier != nullptr) == (work_function != nullptr)) {
            return std::unexpected(at(s.line, "[contact " + s.name + "] needs exactly one of 'barrier' and 'work_function'"));
        }
        c.schottky.by_work_function = work_function != nullptr;
        c.schottky.boundary = patch ? patch->value : std::string{};
    }
    if (region_entry) {  // Unit 68 (MAD6): the electrode is the whole metal region
        if (c.kind != device::ContactKind::electrode) {
            return std::unexpected(at(region_entry->line, "'region' is for a kind = electrode contact"));
        }
        if (patch || !boxes.empty()) {
            return std::unexpected(at(region_entry->line, "an electrode on a metal region takes no 'patch' or 'box'"));
        }
        if (work_function || electrode || metal) {
            return std::unexpected(at(region_entry->line, "the metal region gives the electrode's work function; "
                                                          "drop 'work_function', 'electrode' and 'metal' here"));
        }
        const auto spec = std::ranges::find_if(regions_, [&](const RegionSpec& g) { return g.section->name == region_entry->value; });
        if (spec == regions_.end()) return std::unexpected(at(region_entry->line, "no [region " + region_entry->value + "]"));
        if (!spec->metal) return std::unexpected(at(region_entry->line, "[region " + region_entry->value + "] is not a metal region"));
        c.electrode = *spec->metal;
        const auto r = static_cast<device::RegionId>(spec - regions_.begin());
        for (std::size_t k = 0; k < d.node_region.size(); ++k) {
            if (d.node_region[k] == r) c.nodes.push_back(static_cast<mesh::NodeId>(k));
        }
        if (c.nodes.empty()) return std::unexpected(at(region_entry->line, "metal region '" + region_entry->value + "' has no mesh node"));
        return c;
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
    if (silicide) {  // MAD14: a silicide barrier is measured on silicon
        for (const mesh::NodeId k : c.nodes) {
            const auto& m = d.regions[static_cast<std::size_t>(d.node_region[static_cast<std::size_t>(k)])].material;
            const auto& sp = std::get<physics::Semiconductor>(m).parameters();
            if (!(sp == physics::silicon_parameters) && !(sp == physics::legacy_silicon_parameters)) {
                return std::unexpected(at(silicide->line, "'" + silicide->value + "' is a barrier measured on silicon; "
                                                          "this contact has a node of another material"));
            }
        }
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
    DoseKeys dose;
    for (const Entry& e : s.entries) {
        if (dose_key(e, dose)) {
        } else if (e.key == "transport") {
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
    NITCAD_DECK_TRY(q, dose_charges(dose, std::nullopt));
    if (q.trapped_cm2 != 0.0) f.fixed_charge_cm2 += q.trapped_cm2;  // as for a gate
    if (q.traps_cm2 > 0.0) {
        // Uniform over +-w about E_i: D_it = Delta N_it / (2 w), donor-like below, acceptor-like above.
        const double w = physics::dose_trap_half_width_eV;
        const double dit = q.traps_cm2 / (2.0 * w);
        f.traps.bands.push_back({physics::TrapType::donor, dit, -w, 0.0, 1e-15, 1e-15});
        f.traps.bands.push_back({physics::TrapType::acceptor, dit, 0.0, w, 1e-15, 1e-15});
    }
    return f;
}

bool Builder::dose_key(const Entry& e, DoseKeys& k) {
    const Entry** slot = e.key == "dose"                   ? &k.dose
                         : e.key == "dose_field"           ? &k.field
                         : e.key == "hole_trapping"        ? &k.trapping
                         : e.key == "dose_oxide_thickness" ? &k.thickness
                         : e.key == "interface_trap_yield" ? &k.trap_yield
                                                           : nullptr;
    if (slot) *slot = &e;
    return slot != nullptr;
}

std::expected<Builder::DoseCharges, base::Error> Builder::dose_charges(const DoseKeys& k,
                                                                       std::optional<double> gate_thickness_cm) const {
    if (!k.dose) {
        for (const Entry* e : {k.field, k.trapping, k.thickness, k.trap_yield}) {
            if (e) return std::unexpected(at(e->line, "'" + e->key + "' needs 'dose'"));
        }
        return DoseCharges{};
    }
    NITCAD_DECK_TRY(dose, num(*k.dose));
    if (!std::isfinite(dose) || dose <= 0.0) return std::unexpected(at(k.dose->line, "dose must be finite and positive"));
    if (!k.trapping) return std::unexpected(at(k.dose->line, "dose needs 'hole_trapping'"));
    NITCAD_DECK_TRY(trapping, num(*k.trapping));
    if (!(trapping >= 0.0 && trapping <= 1.0)) {
        return std::unexpected(at(k.trapping->line, "hole_trapping must lie in 0..1"));
    }
    double field = 0.0;
    if (k.field) {
        NITCAD_DECK_TRY(v, num(*k.field));
        if (!std::isfinite(v) || v < 0.0) return std::unexpected(at(k.field->line, "dose_field must be finite and >= 0"));
        field = v;
    }
    double thickness = gate_thickness_cm.value_or(0.0);
    if (!gate_thickness_cm) {
        if (!k.thickness) return std::unexpected(at(k.dose->line, "dose needs 'dose_oxide_thickness'"));
        NITCAD_DECK_TRY(v, len(*k.thickness));
        thickness = v;
    }
    if (!std::isfinite(thickness) || thickness <= 0.0) {
        return std::unexpected(at(k.thickness ? k.thickness->line : k.dose->line,
                                  k.thickness ? "dose_oxide_thickness must be finite and positive"
                                              : "dose needs the gate's oxide_thickness"));
    }
    double yield = 0.0;
    if (k.trap_yield) {
        NITCAD_DECK_TRY(v, num(*k.trap_yield));
        if (!std::isfinite(v) || v < 0.0) {
            return std::unexpected(at(k.trap_yield->line, "interface_trap_yield must be finite and >= 0"));
        }
        yield = v;
    }
    return DoseCharges{physics::oxide_trapped_charge_cm2(dose, field, thickness, trapping), yield * dose};
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

// Unit 68 slice 2 (MAD9; DECISIONS.md amendment 6): with [models] polarization = on, the nitride
// regions' bound polarization charge, lumped into the control volumes. Each region's P (Ga-face,
// along the c-axis, `c_axis = [-]x|y|z`, default y: substrate to surface) is the spontaneous plus
// the piezoelectric part of its strain: the thickest nitride region along c is the relaxed
// reference, every other is pseudomorphic on it, or partly relaxed by its `relaxation = r`
// (0 <= r <= 1; a = a0 + (1 - r)(a_ref - a0)). On each mesh edge i -> j the sheet charge
// (P_i - P_j)(c.e_ij) times the dual face's area goes half to each end, or all to the
// semiconductor end next to an insulator or metal (P 0 there); on the free faces normal to c, P.n
// times the face area (polarization_surface = off: none there, surface donors compensating them).
// polarization_scale = s (0..1) scales all of it.
std::expected<void, base::Error> Builder::polarization(device::DeviceDescription& d) const {
    const Entry *on = nullptr, *scale = nullptr, *axis = nullptr, *surface = nullptr;
    if (const Section* s = single("models")) {
        for (const Entry& e : s->entries) {
            if (e.key == "polarization") on = &e;
            if (e.key == "polarization_scale") scale = &e;
            if (e.key == "c_axis") axis = &e;
            if (e.key == "polarization_surface") surface = &e;
        }
    }
    bool enabled = false;
    if (on) {
        NITCAD_DECK_TRY(v, flag(*on));
        enabled = v;
    }
    for (const Entry* e : {scale, axis, surface}) {
        if (e && !enabled) return std::unexpected(at(e->line, "'" + e->key + "' needs polarization = on"));
    }
    for (const RegionSpec& spec : regions_) {
        if (spec.relaxation && !enabled) {
            return std::unexpected(at(spec.relaxation->line, "'relaxation' needs [models] polarization = on"));
        }
    }
    if (!enabled) return {};
    double k = 1.0;
    if (scale) {
        NITCAD_DECK_TRY(v, num(*scale));
        if (!(v >= 0.0 && v <= 1.0)) return std::unexpected(at(scale->line, "polarization_scale must lie in 0..1"));
        k = v;
    }
    const int dim = d.mesh.dimension();
    std::size_t c = dim == 1 ? 0 : 1;
    double sign = 1.0;
    if (axis) {
        std::string_view v = axis->value;
        if (v.starts_with('-')) {
            sign = -1.0;
            v.remove_prefix(1);
        }
        if (v.size() != 1 || v[0] < 'x' || v[0] > 'z' || v[0] - 'x' >= dim) {
            return std::unexpected(at(axis->line, std::format("c_axis is x, y or z (with '-' to reverse), within this {}D mesh", dim)));
        }
        c = static_cast<std::size_t>(v[0] - 'x');
    }
    // Each region's lattice; a graded region's layers by their own composition.
    const std::size_t R = d.regions.size(), declared = regions_.size();
    std::vector<std::optional<physics::NitrideLattice>> lattice(R);
    std::vector<double> relaxation(R, 0.0);
    for (std::size_t r = 0; r < R; ++r) {
        const RegionSpec& spec = regions_[r < declared ? r : graded_parent_[r - declared]];
        const results::MaterialRecord* record =
            r < declared ? (spec.record ? &*spec.record : nullptr) : &graded_records_[r - declared];
        if (record) {
            std::string_view name = record->material;
            if (const auto user = user_materials_.find(name); user != user_materials_.end()) name = user->second.base;
            lattice[r] = physics::nitride_lattice(name, record->composition.empty() ? 0.0 : record->composition[0]);
        }
        if (spec.relaxation) {
            if (!lattice[r]) return std::unexpected(at(spec.relaxation->line, "'relaxation' is for a nitride region"));
            NITCAD_DECK_TRY(v, num(*spec.relaxation));
            if (!(v >= 0.0 && v <= 1.0)) return std::unexpected(at(spec.relaxation->line, "relaxation must lie in 0..1"));
            relaxation[r] = v;
        }
    }
    // The reference: the thickest nitride region along c (the first of equals).
    const auto points = d.mesh.points();
    std::vector<double> lo(R, std::numeric_limits<double>::infinity()), hi(R, -std::numeric_limits<double>::infinity());
    for (std::size_t i = 0; i < d.node_region.size(); ++i) {
        const auto r = static_cast<std::size_t>(d.node_region[i]);
        lo[r] = std::min(lo[r], points[i][c]);
        hi[r] = std::max(hi[r], points[i][c]);
    }
    std::optional<std::size_t> reference;
    for (std::size_t r = 0; r < R; ++r) {
        if (lattice[r] && (!reference || hi[r] - lo[r] > hi[*reference] - lo[*reference])) reference = r;
    }
    if (!reference) return std::unexpected(at(on->line, "polarization = on, but no region is a nitride"));
    const double a_ref = lattice[*reference]->a0_A;
    std::vector<double> P(R, 0.0);  // [cm^-2], along +c
    for (std::size_t r = 0; r < R; ++r) {
        if (!lattice[r]) continue;
        const double a0 = lattice[r]->a0_A;
        const double a = r == *reference ? a0 : a0 + (1.0 - relaxation[r]) * (a_ref - a0);
        P[r] = k * physics::polarization_sheet_density_cm2(physics::nitride_polarization(*lattice[r], a));
    }
    const auto insulator = [&d](std::size_t node) {
        return device::is_insulator(d.regions[static_cast<std::size_t>(d.node_region[node])]);
    };
    std::vector<double> bound(d.node_region.size(), 0.0);  // charge per node [cm^-2 cm^(D-1)]
    for (const mesh::Edge& e : d.mesh.edges()) {
        const auto i = static_cast<std::size_t>(e.first), j = static_cast<std::size_t>(e.second);
        const double Pi = P[static_cast<std::size_t>(d.node_region[i])], Pj = P[static_cast<std::size_t>(d.node_region[j])];
        if (Pi == Pj) continue;
        const double q = (Pi - Pj) * sign * (points[j][c] - points[i][c]) / e.length * e.coupling_area;
        if (insulator(i)) bound[j] += q;
        else if (insulator(j)) bound[i] += q;
        else bound[i] += 0.5 * q, bound[j] += 0.5 * q;
    }
    bool free_faces = true;  // off: surface donors compensate the free c-faces
    if (surface) {
        NITCAD_DECK_TRY(v, flag(*surface));
        free_faces = v;
    }
    const std::string face = std::string(1, static_cast<char>('x' + c));
    for (const mesh::BoundaryPatch& p : d.mesh.boundary()) {
        const double normal = !free_faces ? 0.0 : p.name == face + "_max" ? sign : p.name == face + "_min" ? -sign : 0.0;
        if (normal == 0.0) continue;
        for (std::size_t n = 0; n < p.nodes.size(); ++n) {
            const auto i = static_cast<std::size_t>(p.nodes[n]);
            bound[i] += P[static_cast<std::size_t>(d.node_region[i])] * normal * p.areas[n];
        }
    }
    const auto volumes = d.mesh.volumes();
    for (std::size_t i = 0; i < bound.size(); ++i) bound[i] /= volumes[i];
    d.bound_charge = std::move(bound);
    return {};
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
        } else if (e.key == "trap_assisted_tunnelling") {
            if (e.value == "off") m.trap_assisted_tunnelling = assemble::TrapTunnelling::off;
            else if (e.value == "hurkx") m.trap_assisted_tunnelling = assemble::TrapTunnelling::hurkx;
            else return std::unexpected(at(e.line, "trap_assisted_tunnelling is off or hurkx"));
        } else if (e.key == "tat_mass_n") {
            NITCAD_DECK_TRY(v, num(e));
            m.tat_mass_n = v;
        } else if (e.key == "tat_mass_p") {
            NITCAD_DECK_TRY(v, num(e));
            m.tat_mass_p = v;
        } else if (e.key == "btbt_nonlocal") {
            if (e.value == "off") m.btbt_nonlocal = assemble::NonlocalTunnelling::off;
            else if (e.value == "kane") m.btbt_nonlocal = assemble::NonlocalTunnelling::kane;
            else if (e.value == "direct_wkb") m.btbt_nonlocal = assemble::NonlocalTunnelling::direct_wkb;
            else return std::unexpected(at(e.line, "btbt_nonlocal is off, kane or direct_wkb"));
        } else if (e.key == "hole_mobility") {
            if (e.value == "caughey_thomas") m.hole_mobility = assemble::HoleMobility::caughey_thomas;
            else if (e.value == "masetti") m.hole_mobility = assemble::HoleMobility::masetti;
            else if (e.value == "automatic") m.hole_mobility = assemble::HoleMobility::automatic;
            else if (e.value == "klaassen") m.hole_mobility = assemble::HoleMobility::klaassen;
            else if (e.value == "klaassen_equilibrium") m.hole_mobility = assemble::HoleMobility::klaassen_equilibrium;
            else return std::unexpected(at(e.line, "hole_mobility is automatic, caughey_thomas, masetti, klaassen or klaassen_equilibrium"));
        } else if (e.key == "electron_mobility") {
            if (e.value == "caughey_thomas") m.electron_mobility = assemble::ElectronMobility::caughey_thomas;
            else if (e.value == "masetti") m.electron_mobility = assemble::ElectronMobility::masetti;
            else if (e.value == "klaassen") m.electron_mobility = assemble::ElectronMobility::klaassen;
            else if (e.value == "klaassen_equilibrium") m.electron_mobility = assemble::ElectronMobility::klaassen_equilibrium;
            else return std::unexpected(at(e.line, "electron_mobility is caughey_thomas, masetti, klaassen or klaassen_equilibrium"));
        } else if (e.key == "surface_mobility") {
            if (e.value == "off") m.surface_mobility = assemble::SurfaceMobility::off;
            else if (e.value == "lombardi") m.surface_mobility = assemble::SurfaceMobility::lombardi;
            else return std::unexpected(at(e.line, "surface_mobility is off or lombardi"));
        } else if (e.key == "quantum_confinement") {
            if (e.value == "off") m.quantum_confinement = assemble::QuantumConfinement::off;
            else if (e.value == "mlda") m.quantum_confinement = assemble::QuantumConfinement::mlda;
            else if (e.value == "density_gradient")
                m.quantum_confinement = assemble::QuantumConfinement::density_gradient;
            else return std::unexpected(at(e.line, "quantum_confinement is off, mlda or density_gradient"));
        } else if (e.key == "mlda_mass_n") {
            NITCAD_DECK_TRY(v, num(e));
            m.mlda_mass_n = v;
        } else if (e.key == "mlda_mass_p") {
            NITCAD_DECK_TRY(v, num(e));
            m.mlda_mass_p = v;
        } else if (e.key == "dg_mass_n") {
            NITCAD_DECK_TRY(v, num(e));
            m.dg_mass_n = v;
        } else if (e.key == "dg_mass_p") {
            NITCAD_DECK_TRY(v, num(e));
            m.dg_mass_p = v;
        } else if (e.key == "dg_gamma_n") {
            NITCAD_DECK_TRY(v, num(e));
            m.dg_gamma_n = v;
        } else if (e.key == "dg_gamma_p") {
            NITCAD_DECK_TRY(v, num(e));
            m.dg_gamma_p = v;
        } else if (e.key == "ionization_levels") {
            if (e.value == "fixed") m.ionization_levels = assemble::IonizationLevels::fixed;
            else if (e.value == "altermatt") m.ionization_levels = assemble::IonizationLevels::altermatt;
            else return std::unexpected(at(e.line, "ionization_levels is fixed or altermatt"));
        } else if (e.key == "impact_current_resolution") {
            NITCAD_DECK_TRY(v, num(e));
            m.impact_current_resolution = v;
        } else if (e.key == "polarization" || e.key == "polarization_scale" || e.key == "c_axis" ||
                   e.key == "polarization_surface") {
            // Unit 68 slice 2: read with the device (Builder::polarization).
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
        } else if (e.key == "natural_norm") {  // Unit 59 (NQ7)
            if (e.value == "convergence") {
                o.newton.line_search.norm = solve::NaturalNorm::convergence;
            } else if (e.value == "weighted") {
                o.newton.line_search.norm = solve::NaturalNorm::weighted;
            } else {
                return std::unexpected(at(e.line, "natural_norm is convergence or weighted"));
            }
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
            // Unit 58 SQ7 (a): gmres is gmres_ilut, its preconditioner set apart (an unchanged identity).
            else if (e.value == "gmres_ilut" || e.value == "gmres")
                o.linear.backend = linalg::SolverBackend::gmres_ilut;
            else return std::unexpected(at(e.line, "backend is eigen, pardiso, gmres or gmres_ilut"));
        } else if (e.key == "preconditioner") {  // Unit 58 SD6
            if (e.value == "ilut") o.linear.preconditioner = linalg::Preconditioner::ilut;
            else if (e.value == "amg") o.linear.preconditioner = linalg::Preconditioner::amg;
            else if (e.value == "cpr") o.linear.preconditioner = linalg::Preconditioner::cpr;
            else if (e.value == "cpr_amg") o.linear.preconditioner = linalg::Preconditioner::cpr_amg;  // Unit 63 MD4
            else return std::unexpected(at(e.line, "preconditioner is ilut, amg, cpr or cpr_amg"));
        } else if (e.key == "amg_strength_threshold") {
            NITCAD_DECK_TRY(v, num(e));
            if (!(v >= 0.0 && v < 1.0)) return std::unexpected(at(e.line, "amg_strength_threshold is in [0, 1)"));
            o.linear.amg_strength_threshold = v;
        } else if (e.key == "amg_smoother_degree") {
            NITCAD_DECK_TRY(v, count(e, 1, 10));
            o.linear.amg_smoother_degree = v;
        } else if (e.key == "amg_coarse_size") {
            NITCAD_DECK_TRY(v, count(e, 1, 100000));
            o.linear.amg_coarse_size = v;
        } else if (e.key == "amg_max_levels") {
            NITCAD_DECK_TRY(v, count(e, 2, 50));
            o.linear.amg_max_levels = v;
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
    if (o.linear.preconditioner != linalg::Preconditioner::ilut &&
        o.linear.backend != linalg::SolverBackend::gmres_ilut) {
        return std::unexpected(at(s->line, "preconditioner = amg, cpr or cpr_amg needs backend = gmres"));
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
    const Entry* initial = nullptr;
    const Entry* restart = nullptr;
    NITCAD_DECK_CHECK(operating_points(s, plan, [&](const Entry& e) -> std::expected<bool, base::Error> {
        if (e.key == "restart") restart = &e;
        else if (e.key == "initial") initial = &e;
        else return false;
        return true;
    }));
    return starts(initial, restart, plan);
}

// `initial` and `restart` of an analysis (Unit 61, RD2): the files, relative to the deck. Both in one
// section is an error at the later one's line.
std::expected<void, base::Error> Builder::starts(const Entry* initial, const Entry* restart, Plan& plan) const {
    if (initial != nullptr && restart != nullptr) {
        return std::unexpected(at(std::max(initial->line, restart->line),
                                  "initial and restart together: restart continues a run, initial only guesses its "
                                  "first state"));
    }
    if (initial != nullptr) plan.initial = file(*initial);
    if (restart != nullptr) plan.restart = file(*restart);
    return {};
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
    const Entry* initial = nullptr;
    const Entry* restart = nullptr;
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
            NITCAD_DECK_TRY(v, lift(e.line, seconds(e.value)));
            o.t_end_s = v;
            t_end = true;
        } else if (e.key == "dt_initial" || e.key == "dt_min" || e.key == "dt_max") {
            NITCAD_DECK_TRY(v, lift(e.line, seconds(e.value)));
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
            NITCAD_DECK_TRY(v, lift(e.line, series(e.value, [](std::string_view t) { return seconds(t); })));
            o.output_times_s = std::move(v);
        } else if (e.key == "initial") {
            initial = &e;
        } else if (e.key == "restart") {
            restart = &e;
        } else {
            return std::unexpected(unknown(e, s));
        }
    }
    if (!t_end) return std::unexpected(at(s.line, "[transient] needs t_end"));
    NITCAD_DECK_CHECK(starts(initial, restart, plan));
    for (auto& f : waveforms) plan.waveforms.push_back(f ? std::move(*f) : solve::Waveform::constant(0.0));
    return {};
}

// [heavy_ion <name>] sections (Unit 65, HD9): one track each, a [transient] option only. Lengths as
// every deck length (cm, or with um or nm); `area` in cm^2 or with um2; `let` in pC/um, `let_mev`
// in MeV cm^2/mg (exactly one); times with their units. No defaults (HQ4): start, direction,
// length, let or let_mev, t0 and sigma always; radius in 2D and 3D; depth in 2D; area in 1D.
std::expected<void, base::Error> Builder::heavy_ions(Plan& plan) const {
    const int dim = plan.device.mesh().dimension();
    for (const Section& s : d_.sections) {
        if (s.kind != "heavy_ion") continue;
        if (plan.analysis != Analysis::transient) {
            return std::unexpected(at(s.line, "[heavy_ion " + s.name + "] needs a [transient] analysis"));
        }
        physics::HeavyIonTrack t;
        std::set<std::string> given;
        for (const Entry& e : s.entries) {
            if (!given.insert(e.key).second) return std::unexpected(at(e.line, "a second '" + e.key + "'"));
            if (e.key == "start") {
                NITCAD_DECK_TRY(v, lengths(e, e.value, static_cast<std::size_t>(dim)));
                for (std::size_t k = 0; k < v.size(); ++k) t.start_cm[k] = v[k];
            } else if (e.key == "direction") {
                std::vector<double> v;
                for (const std::string_view item : items(e.value)) {
                    for (const std::string_view w : tokens(item)) {
                        NITCAD_DECK_TRY(x, lift(e.line, number(w)));
                        v.push_back(x);
                    }
                }
                if (v.size() != static_cast<std::size_t>(dim)) {
                    return std::unexpected(at(e.line, "'direction' needs " + std::to_string(dim) + " components"));
                }
                for (std::size_t k = 0; k < v.size(); ++k) t.direction[k] = v[k];
            } else if (e.key == "length" || e.key == "radius" || e.key == "depth") {
                NITCAD_DECK_TRY(v, len(e));
                (e.key == "length" ? t.length_cm : e.key == "radius" ? t.radius_cm : t.depth_cm) = v;
            } else if (e.key == "area") {
                std::string_view text = trimmed(e.value);
                double scale = 1.0;
                if (text.ends_with("um2")) {
                    text.remove_suffix(3);
                    scale = 1e-8;
                }
                NITCAD_DECK_TRY(v, lift(e.line, number(text)));
                t.area_cm2 = v * scale;
            } else if (e.key == "let") {
                NITCAD_DECK_TRY(v, num(e));
                t.let_C_per_cm = v * 1e-12 / 1e-4;  // pC/um to C/cm
            } else if (e.key == "let_mev") {
                NITCAD_DECK_TRY(v, num(e));
                t.let_C_per_cm = physics::let_from_mev_cm2_mg(v);
            } else if (e.key == "t0" || e.key == "sigma") {
                NITCAD_DECK_TRY(v, lift(e.line, seconds(e.value)));
                (e.key == "t0" ? t.t0_s : t.sigma_s) = v;
            } else {
                return std::unexpected(unknown(e, s));
            }
        }
        if (given.contains("let") && given.contains("let_mev")) {
            return std::unexpected(at(s.line, "[heavy_ion " + s.name + "] takes 'let' or 'let_mev', not both"));
        }
        std::vector<std::string> needed{"start", "direction", "length", "t0", "sigma"};
        if (dim >= 2) needed.emplace_back("radius");
        if (dim == 2) needed.emplace_back("depth");
        if (dim == 1) needed.emplace_back("area");
        for (const std::string& k : needed) {
            if (!given.contains(k)) return std::unexpected(at(s.line, "[heavy_ion " + s.name + "] needs '" + k + "'"));
        }
        if (!given.contains("let") && !given.contains("let_mev")) {
            return std::unexpected(at(s.line, "[heavy_ion " + s.name + "] needs 'let' (pC/um) or 'let_mev'"));
        }
        if (auto ok = physics::check_heavy_ion_track(t, dim); !ok) {
            return std::unexpected(at(s.line, "[heavy_ion " + s.name + "]: " + ok.error().message));
        }
        plan.transient.heavy_ions.push_back(t);
    }
    return {};
}

std::expected<void, base::Error> Builder::small_signal(const Section& s, Plan& plan) const {
    plan.analysis = Analysis::small_signal;
    bool frequencies = false;
    const Entry* initial = nullptr;
    const Entry* restart = nullptr;
    NITCAD_DECK_CHECK(operating_points(s, plan, [&](const Entry& e) -> std::expected<bool, base::Error> {
        if (e.key == "restart") {
            restart = &e;
            return true;
        }
        if (e.key == "initial") {
            initial = &e;
            return true;
        }
        if (e.key != "frequencies" && e.key != "field_frequencies") return false;
        NITCAD_DECK_TRY(v, lift(e.line, series(e.value, [](std::string_view t) { return hertz(t); })));
        if (e.key == "frequencies") {
            plan.small_signal.frequencies_Hz = std::move(v);
            frequencies = true;
        } else {
            plan.small_signal.field_frequencies_Hz = std::move(v);
        }
        return true;
    }));
    if (!frequencies) return std::unexpected(at(s.line, "[small_signal] needs frequencies"));
    NITCAD_DECK_CHECK(starts(initial, restart, plan));
    return {};
}

std::expected<void, base::Error> Builder::trace(const Section& s, Plan& plan) const {
    plan.analysis = Analysis::trace;
    solve::TraceOptions& o = plan.trace;
    o.start_V.assign(plan.device.contacts().size(), 0.0);
    std::vector<bool> given(o.start_V.size(), false), given_t;
    bool has_contact = false, has_end = false;
    const Entry* initial = nullptr;
    const Entry* restart = nullptr;
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
        } else if (e.key == "initial") {
            initial = &e;
        } else if (e.key == "restart") {
            restart = &e;
        } else {
            return std::unexpected(unknown(e, s));
        }
    }
    if (!has_contact) return std::unexpected(at(s.line, "[trace] needs contact, the swept contact"));
    if (!has_end) return std::unexpected(at(s.line, "[trace] needs end, the swept contact's last bias"));
    NITCAD_DECK_CHECK(starts(initial, restart, plan));
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
        p.reserve(bias.size());
        for (const auto& v : bias) p.push_back(at_point(v, k));
        plan.points.push_back(std::move(p));
        if (any_thermal) {
            std::vector<double> t;
            t.reserve(temperature.size());
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
                                        : e.key == "adapt"        ? &plan.outputs.adapt
                                        : e.key == "adapt_meshes" ? &plan.outputs.adapt_meshes
                                        : e.key == "adapt_history" ? &plan.outputs.adapt_history
                                                                  : nullptr;
        if (!target) return std::unexpected(unknown(e, *s));
        if ((e.key == "adapt" || e.key == "adapt_meshes" || e.key == "adapt_history") && !plan.adapt) {
            return std::unexpected(at(e.line, "'" + e.key + "' is written by an [adapt] run; the deck has no [adapt]"));
        }
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
    return {};
}

// [adapt] (Unit 17, AD7; DECISIONS.md AR8): the loop's options, checked with the line of the key that
// is wrong; the loop's own checks (solve::adapt) repeat these without lines.
std::expected<void, base::Error> Builder::adapt(Plan& plan) const {
    const Section* s = single("adapt");
    if (!s) return {};
    if (plan.analysis != Analysis::sweep) {
        return std::unexpected(at(s->line, "[adapt] refines the mesh of a [sweep] (one mesh for the sweep); this deck's "
                                           "analysis is not a sweep"));
    }
    if (!plan.restart.empty() || !plan.initial.empty()) {
        return std::unexpected(at(s->line, "[adapt] with a restart or initial: the checkpoint belongs to the deck's mesh, not the "
                                           "adapted one"));
    }
    if (type_ == MeshType::generate3d || type_ == MeshType::file || (type_ == MeshType::extrude3d && !swept_)) {
        return std::unexpected(at(s->line, "[adapt] refines a tensor, generate2d or swept extrude3d (sweep = on) mesh; "
                                           "a general 3D or file mesh is not adapted in this unit (AQ2)"));
    }
    AdaptPlan a = mesh_info_;
    if (type_ == MeshType::tensor) a.kind = solve::AdaptMeshKind::tensor;
    solve::AdaptOptions& o = a.options;
    const Entry* weighting = nullptr;
    const Entry* mode = nullptr;
    const Entry* budget_dofs = nullptr;
    const Entry* budget_seconds = nullptr;

    const Entry* max_points = nullptr;
    const Entry* predictor = nullptr;
    // Unit 17c: the first derived goal, an explicit `at`, `vt_tolerance` (for their lines).
    const Entry* derived_goal = nullptr;
    const Entry* at_entry = nullptr;
    const Entry* vt_tolerance = nullptr;
    bool has_vt = false;
    std::set<std::string> seen;
    for (const Entry& e : s->entries) {
        if (e.key != "goal" && !seen.insert(e.key).second) return std::unexpected(at(e.line, "a second '" + e.key + "'"));
        if (e.key == "goal") {
            const auto w = tokens(e.value);
            // Unit 17c (DG1): vt, gm and swing of the contact's transfer curve beside its current.
            solve::GoalKind goal_kind = solve::GoalKind::current;
            const bool known = !w.empty() && (w[0] == "current" || w[0] == "vt" || w[0] == "gm" || w[0] == "swing");
            if (!known || w.size() != 2) {
                return std::unexpected(at(e.line, "unknown goal '" + e.value +
                                                      "': a goal is 'current <contact>', 'vt <contact>', 'gm <contact>' "
                                                      "or 'swing <contact>'"));
            }
            if (w[0] == "vt") goal_kind = solve::GoalKind::vt;
            else if (w[0] == "gm") goal_kind = solve::GoalKind::gm;
            else if (w[0] == "swing") goal_kind = solve::GoalKind::swing;
            NITCAD_DECK_TRY(c, contact_index(e, w[1], plan));
            const auto kind = plan.device.contacts()[c].kind;
            if (kind == device::ContactKind::gate || kind == device::ContactKind::electrode) {
                return std::unexpected(at(e.line, "the goal's contact '" + std::string(w[1]) +
                                                      "' is a gate or an electrode, which carries no current"));
            }
            for (const solve::AdaptGoal& g : o.goals) {
                if (g.contact == c && g.kind == goal_kind) {
                    // A current goal's message as before Unit 17c.
                    const std::string kind_word = goal_kind == solve::GoalKind::current ? "" : std::string(w[0]) + " ";
                    return std::unexpected(at(e.line, "a second " + kind_word + "goal on '" + std::string(w[1]) + "'"));
                }
            }
            if (goal_kind != solve::GoalKind::current && derived_goal == nullptr) derived_goal = &e;
            if (goal_kind == solve::GoalKind::vt) has_vt = true;
            o.goals.push_back({.contact = c, .kind = goal_kind});
        } else if (e.key == "tolerance") {
            NITCAD_DECK_TRY(v, num(e));
            if (!(v > 0.0 && v < 1.0)) return std::unexpected(at(e.line, "the tolerance is relative, in (0, 1)"));
            o.tolerance = v;
        } else if (e.key == "vt_tolerance") {
            // DG5: a vt goal's absolute allowance, in volts.
            vt_tolerance = &e;
            NITCAD_DECK_TRY(v, num(e));
            if (!(v > 0.0 && std::isfinite(v))) return std::unexpected(at(e.line, "vt_tolerance must be a positive voltage"));
            o.vt_tolerance = v;
        } else if (e.key == "at") {
            at_entry = &e;
            if (e.value == "auto") continue;
            for (const std::string_view item : items(e.value)) {
                NITCAD_DECK_TRY(k, lift(e.line, integer(item)));
                if (k < 0 || static_cast<std::size_t>(k) >= plan.points.size()) {
                    return std::unexpected(at(e.line, "the adaptation point " + std::to_string(k) + " is not a point of the sweep (0.." +
                                                          std::to_string(plan.points.size() - 1) + ")"));
                }
                o.at.push_back(static_cast<std::size_t>(k));
            }
        } else if (e.key == "max_points") {
            max_points = &e;
            NITCAD_DECK_TRY(v, count(e, 1, 1000));
            o.intelligent.max_points = static_cast<std::size_t>(v);
        } else if (e.key == "indicators") {
            if (e.value == "auto") continue;
            o.indicators.enabled.fill(false);
            for (const std::string_view item : items(e.value)) {
                std::size_t k = 0;
                while (k < solve::indicator_count && item != solve::to_string(static_cast<solve::Indicator>(k))) ++k;
                if (k == solve::indicator_count) {
                    return std::unexpected(at(e.line, "unknown indicator '" + std::string(item) +
                                                          "' (potential, electron_density, hole_density, doping, generation)"));
                }
                o.indicators.enabled[k] = true;
            }
        } else if (e.key == "weighting") {
            weighting = &e;
            if (e.value == "adjoint") o.weighting = solve::Weighting::adjoint;
            else if (e.value != "none") return std::unexpected(at(e.line, "weighting is none or adjoint"));
        } else if (e.key == "policy") {
            if (e.value == "intelligent") o.policy = solve::Policy::intelligent;
            else if (e.value != "basic") return std::unexpected(at(e.line, "policy is basic or intelligent"));
        } else if (e.key == "mode") {
            if (e.value != "tolerance" && e.value != "budget") return std::unexpected(at(e.line, "mode is tolerance or budget"));
            mode = &e;
        } else if (e.key == "budget_dofs") {
            budget_dofs = &e;
            NITCAD_DECK_TRY(v, count(e, 1, std::numeric_limits<int>::max()));
            o.intelligent.dof_budget = static_cast<std::size_t>(v);
        } else if (e.key == "budget_seconds") {
            budget_seconds = &e;
            NITCAD_DECK_TRY(v, num(e));
            if (!(v > 0.0 && std::isfinite(v))) return std::unexpected(at(e.line, "budget_seconds must be positive"));
            o.intelligent.seconds_budget = v;
        } else if (e.key == "verify") {
            if (e.value == "off") o.verify = false;
            else if (e.value != "richardson") return std::unexpected(at(e.line, "verify is richardson or off"));
        } else if (e.key == "certify") {
            if (e.value == "required") a.certify_required = true;
            else if (e.value != "report") return std::unexpected(at(e.line, "certify is report or required"));
        } else if (e.key == "max_passes") {
            NITCAD_DECK_TRY(v, count(e, 1, 1000));
            o.max_passes = v;
        } else if (e.key == "max_nodes") {
            NITCAD_DECK_TRY(v, count(e, 1, std::numeric_limits<int>::max()));
            o.max_nodes = static_cast<std::size_t>(v);
        } else if (e.key == "min_size") {
            NITCAD_DECK_TRY(v, len(e));
            if (!(v > 0.0 && std::isfinite(v))) return std::unexpected(at(e.line, "min_size must be a positive length"));
            o.min_size = v;
        } else if (e.key == "grading") {
            NITCAD_DECK_TRY(v, num(e));
            if (!(v >= 1.0 && std::isfinite(v))) return std::unexpected(at(e.line, "grading must be >= 1"));
            o.grading = v;
        } else if (e.key == "theta") {
            NITCAD_DECK_TRY(v, num(e));
            if (!(v > 0.0 && v <= 1.0)) return std::unexpected(at(e.line, "theta must be in (0, 1]"));
            o.theta = v;
        } else if (e.key == "predictor") {
            predictor = &e;
        } else {
            return std::unexpected(unknown(e, *s));
        }
    }
    if (o.goals.empty()) return std::unexpected(at(s->line, "[adapt] needs at least one 'goal = current <contact>'"));
    // Unit 17c (DG1, DG2, DG5): a derived goal's sweep, points and tolerance.
    if (vt_tolerance != nullptr && !has_vt) {
        return std::unexpected(at(vt_tolerance->line, "vt_tolerance is for a 'goal = vt <contact>'"));
    }
    if (derived_goal != nullptr) {
        if (at_entry != nullptr && !o.at.empty()) {
            return std::unexpected(at(at_entry->line, "adaptation points cannot be given with a vt, gm or swing goal "
                                                      "(its window on the transfer curve decides; use at = auto)"));
        }
        if (plan.points.size() < 3) {
            return std::unexpected(at(derived_goal->line, "a vt, gm or swing goal needs a sweep of at least 3 points"));
        }
        const std::vector<std::size_t> swept = solve::swept_contacts(plan.points);  // DR9: adapt()'s rule
        if (swept.size() != 1) {
            return std::unexpected(at(derived_goal->line, "a vt, gm or swing goal needs a sweep of exactly one contact "
                                                          "(this sweep moves " + std::to_string(swept.size()) + ")"));
        }
        const auto& swept_contact = plan.device.contacts()[swept[0]];
        if (swept_contact.kind != device::ContactKind::gate && swept_contact.kind != device::ContactKind::electrode) {
            return std::unexpected(at(derived_goal->line, "a vt, gm or swing goal needs the swept contact to be a gate or "
                                                          "an electrode; '" + swept_contact.name + "' is not"));
        }
    }
    const bool budget = mode != nullptr && mode->value == "budget";
    if (budget && o.policy != solve::Policy::intelligent) {
        return std::unexpected(at(mode->line, "mode = budget needs policy = intelligent (the basic loop takes max_nodes)"));
    }
    if (budget && budget_dofs == nullptr && budget_seconds == nullptr) {
        return std::unexpected(at(mode->line, "mode = budget needs budget_dofs or budget_seconds"));
    }
    if (!budget && (budget_dofs != nullptr || budget_seconds != nullptr)) {
        const Entry* e = budget_dofs != nullptr ? budget_dofs : budget_seconds;
        return std::unexpected(at(e->line, "'" + e->key + "' is for mode = budget"));
    }
    // at = auto (or absent): the first and last point, or with the intelligent policy up to
    // max_points chosen from a coarse sweep (ID3; default 3, AR8).
    if (o.at.empty() && o.policy == solve::Policy::intelligent && o.intelligent.max_points == 0) o.intelligent.max_points = 3;
    if (max_points != nullptr && (o.policy != solve::Policy::intelligent || !o.at.empty())) {
        return std::unexpected(at(max_points->line, "max_points is for at = auto with policy = intelligent"));
    }
    if (o.weighting == solve::Weighting::adjoint) {
        if (auto e = solve::check_adjoint(plan.options)) {
            e->message = d_.source + ":" + std::to_string(weighting->line) + ": weighting = adjoint: " + e->message;
            e->context = base::ErrorContext{.index = static_cast<std::size_t>(weighting->line), .value = std::nullopt};
            return std::unexpected(std::move(*e));
        }
        for (const solve::AdaptGoal& g : o.goals) {
            if (plan.device.contacts()[g.contact].kind != device::ContactKind::ohmic) {
                return std::unexpected(at(weighting->line, "weighting = adjoint needs ohmic goal contacts; '" +
                                                               plan.device.contacts()[g.contact].name + "' is not"));
            }
        }
    }
    // Unit 17b (LB9): the learned model, loaded now so its versions are checked by `nitcad check`.
    if (predictor != nullptr) {
        if (o.policy != solve::Policy::intelligent) {
            return std::unexpected(at(predictor->line, "predictor needs policy = intelligent (it replaces the a priori prediction)"));
        }
        a.predictor_path = file(*predictor);
        auto model = io::load_learned_model(a.predictor_path);
        if (!model) return std::unexpected(at(predictor->line, "predictor: " + model.error().message));
        auto digest = io::file_digest(a.predictor_path);
        if (!digest) return std::unexpected(at(predictor->line, "predictor: " + digest.error().message));
        a.predictor = std::make_shared<const results::LearnedModel>(std::move(*model));
        a.predictor_digest = *digest;
    }
    plan.adapt = std::move(a);
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
    NITCAD_DECK_CHECK(user_materials());
    for (const Section& s : d_.sections) {
        if (s.kind != "region") continue;
        NITCAD_DECK_TRY(spec, region(s));
        regions_.push_back(std::move(spec));
    }
    if (regions_.empty()) return std::unexpected(at(1, "the deck declares no [region]"));
    for (const RegionSpec& g : regions_) d.regions.push_back(*g.region);
    NITCAD_DECK_CHECK(build_mesh(d));
    NITCAD_DECK_CHECK(doping(d));
    NITCAD_DECK_CHECK(grade(d));
    NITCAD_DECK_CHECK(polarization(d));
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
    Plan plan{.device = std::move(*device), .generated = std::move(generated_), .mesh_type = type_name_};
    for (const RegionSpec& g : regions_) {
        if (g.record) plan.materials.push_back(*g.record);
    }
    plan.materials.insert(plan.materials.end(), graded_records_.begin(), graded_records_.end());
    NITCAD_DECK_CHECK(models(plan.options.models));
    NITCAD_DECK_CHECK(solver(plan));
    NITCAD_DECK_CHECK(analysis(plan));
    NITCAD_DECK_CHECK(heavy_ions(plan));
    NITCAD_DECK_CHECK(adapt(plan));
    NITCAD_DECK_CHECK(outputs(plan));
    return plan;
}

}  // namespace

std::expected<device::DopingProfile, std::string> parse_profile(std::string_view value) {
    using namespace detail;
    const auto t = tokens(value);
    if (t.size() < 3) return std::unexpected(std::string("a profile starts '<donor|acceptor> <shape> <peak>'"));
    device::DopingProfile p;
    if (t[0] == "donor") p.dopant = device::Dopant::donor;
    else if (t[0] == "acceptor") p.dopant = device::Dopant::acceptor;
    else return std::unexpected(std::string("the dopant is donor or acceptor"));
    if (t[1] == "constant") p.shape = device::DepthShape::constant;
    else if (t[1] == "gaussian") p.shape = device::DepthShape::gaussian;
    else if (t[1] == "erfc") p.shape = device::DepthShape::erfc;
    else return std::unexpected(std::string("the shape is constant, gaussian or erfc"));
    const auto peak = number(t[2]);
    if (!peak) return std::unexpected(peak.error());
    p.peak_cm3 = *peak;
    for (std::size_t k = 3; k < t.size(); ++k) {
        const std::size_t eq = t[k].find('=');
        if (eq == std::string_view::npos) return std::unexpected("expected name=value, not '" + std::string(t[k]) + "'");
        const std::string_view name = t[k].substr(0, eq);
        const std::string_view v = t[k].substr(eq + 1);
        static constexpr std::array<std::string_view, 3> axes{"x", "y", "z"};
        const auto a = std::find(axes.begin(), axes.end(), name);
        if (name == "axis") {
            const auto w = std::find(axes.begin(), axes.end(), v);
            if (w == axes.end()) return std::unexpected(std::string("axis is x, y or z"));
            p.axis = static_cast<int>(w - axes.begin());
        } else if (name == "center" || name == "width" || name == "lateral") {
            const auto l = length(v);
            if (!l) return std::unexpected(l.error());
            (name == "center" ? p.center : name == "width" ? p.width : p.lateral) = *l;
        } else if (a != axes.end()) {
            const auto range = items(v, ':');
            if (range.size() != 2) return std::unexpected(std::string("a window is <lo>:<hi>"));
            const auto lo = length(range[0]);
            if (!lo) return std::unexpected(lo.error());
            const auto hi = length(range[1]);
            if (!hi) return std::unexpected(hi.error());
            const auto k3 = static_cast<std::size_t>(a - axes.begin());
            p.lo[k3] = *lo;
            p.hi[k3] = *hi;
        } else {
            return std::unexpected("unknown profile attribute '" + std::string(name) + "'");
        }
    }
    return p;
}

std::string profile_text(const device::DopingProfile& p) {
    // Lengths in cm, shortest round-trip: the same doubles read back.
    static constexpr const char* axes[] = {"x", "y", "z"};
    std::string out = std::string(p.dopant == device::Dopant::donor ? "donor " : "acceptor ") +
                      (p.shape == device::DepthShape::constant ? "constant " : p.shape == device::DepthShape::gaussian ? "gaussian " : "erfc ") +
                      std::format("{}", p.peak_cm3);
    if (p.shape != device::DepthShape::constant || p.axis != 0) out += std::string(" axis=") + axes[p.axis];
    if (p.shape != device::DepthShape::constant) out += std::format(" center={}cm width={}cm", p.center, p.width);
    for (int k = 0; k < 3; ++k) {
        if (k == p.axis && p.shape != device::DepthShape::constant) continue;
        const auto i = static_cast<std::size_t>(k);
        if (std::isinf(p.lo[i]) && std::isinf(p.hi[i])) continue;
        out += std::format(" {}={}cm:{}cm", axes[k], p.lo[i], p.hi[i]);
    }
    if (p.lateral != 0.0) out += std::format(" lateral={}cm", p.lateral);
    return out;
}


// Unit 26m MT8: the parameters substituted first, so the build reads values as written out (the plan keeps no view
// into the deck it read; [adapt] keeps the deck as given and substitutes again per pass).
std::expected<Plan, base::Error> build(const Deck& deck) {
    const auto read = substituted(deck);
    if (!read) return std::unexpected(read.error());
    auto plan = Builder(*read, nullptr).run();
    if (plan && plan->adapt) plan->adapt->deck = std::make_shared<const Deck>(deck);
    return plan;
}

std::expected<Plan, base::Error> build(const Deck& deck, const solve::MeshRefinement& refinement) {
    const auto read = substituted(deck);
    if (!read) return std::unexpected(read.error());
    return Builder(*read, &refinement).run();
}

solve::AdaptableMesh adaptable(Plan plan) {
    solve::AdaptableMesh m{.device = std::move(plan.device)};
    if (plan.adapt) {
        const AdaptPlan& a = *plan.adapt;
        m.kind = a.kind;
        m.axis_protected = a.axis_protected;
        m.sweep_axis = a.sweep_axis;
        m.sweep_protected = a.sweep_protected;
        m.interface_resolution = a.interface_resolution;
        m.max_size = a.max_size;
    }
    if (plan.generated) {
        m.section = plan.generated->section;
        m.generated = std::move(plan.generated);
    }
    return m;
}

solve::DeviceBuilder device_builder(const Plan& plan) {
    std::shared_ptr<const Deck> deck = plan.adapt ? plan.adapt->deck : nullptr;
    return [deck](const solve::MeshRefinement& r) -> std::expected<solve::AdaptableMesh, base::Error> {
        if (!deck) return std::unexpected(base::Error{base::ErrorCode::invalid_input, "the plan has no [adapt]", std::nullopt});
        auto p = build(*deck, r);
        if (!p) return std::unexpected(std::move(p.error()));
        return adaptable(std::move(*p));
    };
}

}  // namespace NiTCAD::deck

#undef NITCAD_DECK_TRY
#undef NITCAD_DECK_CHECK
