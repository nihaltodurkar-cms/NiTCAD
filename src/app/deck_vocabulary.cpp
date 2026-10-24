// The deck's vocabulary (deck_vocabulary.hpp), transcribed from deck/build.hpp's reference and the
// keys deck::build reads; the vocabulary test checks both directions against deck::build.
#include "NiTCAD/app/deck_vocabulary.hpp"

#include <algorithm>
#include <array>
#include <numeric>

#include "NiTCAD/app/deck_sections.hpp"

namespace NiTCAD::app {

namespace {

using Words = std::span<const std::string_view>;

constexpr std::array<std::string_view, 2> on_off{"on", "off"};
constexpr std::array<std::string_view, 4> mesh_types{"tensor", "generate2d", "generate3d", "file"};
constexpr std::array<std::string_view, 9> materials{"silicon", "silicon_legacy", "germanium", "gaas", "ingaas",
                                                    "sic_4h",  "algaas",         "sio2",      "insulator"};
constexpr std::array<std::string_view, 2> dopants{"donor", "acceptor"};
constexpr std::array<std::string_view, 3> shapes{"constant", "gaussian", "erfc"};
constexpr std::array<std::string_view, 4> contact_kinds{"ohmic", "gate", "electrode", "schottky"};
constexpr std::array<std::string_view, 3> electrodes{"n_poly", "p_poly", "metal"};
constexpr std::array<std::string_view, 2> transports{"drift_diffusion", "thermionic"};
constexpr std::array<std::string_view, 2> thermal_kinds{"isothermal", "resistance"};
constexpr std::array<std::string_view, 3> btbt_nonlocal{"off", "kane", "direct_wkb"};
constexpr std::array<std::string_view, 3> hole_mobility{"automatic", "caughey_thomas", "masetti"};
constexpr std::array<std::string_view, 2> electron_mobility{"caughey_thomas", "masetti"};
constexpr std::array<std::string_view, 2> surface_mobility{"off", "lombardi"};
constexpr std::array<std::string_view, 2> quantum_confinement{"off", "mlda"};
constexpr std::array<std::string_view, 2> ionization_levels{"fixed", "altermatt"};
constexpr std::array<std::string_view, 2> equations{"drift_diffusion", "equilibrium_poisson"};
constexpr std::array<std::string_view, 3> line_search{"off", "residual", "natural"};
constexpr std::array<std::string_view, 3> backends{"eigen", "pardiso", "gmres_ilut"};
constexpr std::array<std::string_view, 2> fallbacks{"sparse_lu", "none"};
constexpr std::array<std::string_view, 2> integrators{"bdf2", "backward_euler"};
constexpr std::array<std::string_view, 6> waveforms{"constant", "step", "ramp", "pulse", "pwl", "sine"};

constexpr std::array<KeyInfo, 1> device_keys{{{"temperature", 0, false, "300", "lattice temperature, K", {}, {}, {}}}};

constexpr std::array<KeyInfo, 18> mesh_keys{{
    {"type", 0, false, "tensor", "tensor | generate2d | generate3d | file", mesh_types, {}, {}},
    {"x", 0, false, "<axis>", "tensor: x nodes, lengths or linspace(a, b, n) pieces joined by +", {}, {}, "tensor"},
    {"y", 0, false, "<axis>", "tensor: y nodes", {}, {}, "tensor"},
    {"z", 0, false, "<axis>", "tensor: z nodes", {}, {}, "tensor"},
    {"straddle", 1, true, "<position> <spacing>", "tensor: a material interface on axis x, y or z", {}, {}, "tensor"},
    {"patch", 1, true, "<x y>, <x y>", "generate2d: a boundary patch", {}, {}, "generate2d"},
    {"patch", 2, true, "<x y z>, ...", "generate3d: patch <face> <name>, a polygon tiling the face", {}, {},
     "generate3d"},
    {"max_size", 0, false, "<length>", "generated meshes: largest element size", {}, {}, "generate2d generate3d"},
    {"interface_resolution", 0, false, "<length>", "generated meshes: element size at interfaces", {}, {},
     "generate2d generate3d"},
    {"min_angle", 0, false, "20.7", "generate2d: smallest triangle angle, degrees", {}, {}, "generate2d"},
    {"radius_edge", 0, false, "2", "generate3d: largest radius-edge ratio", {}, {}, "generate3d"},
    {"max_nodes", 0, false, "2000000", "generated meshes: node limit", {}, {}, "generate2d generate3d"},
    {"box", 0, true, "<lo> <hi> <h>", "generated meshes: a smaller size h inside a box", {}, {},
     "generate2d generate3d"},
    {"vertex", 0, true, "<x y z>", "generate3d: a domain vertex, numbered from 0", {}, {}, "generate3d"},
    {"face", 1, true, "<i j k ...>", "generate3d: a domain face by its vertices, also a patch", {}, {},
     "generate3d"},
    {"cuboid", 0, false, "<x0 x1 y0 y1 z0 z1>", "generate3d: a box domain with faces x_min ... z_max", {}, {},
     "generate3d"},
    {"plane", 0, true, "<x y z>, <x y z>, <x y z>", "generate3d: a planar interface", {}, {}, "generate3d"},
    {"path", 0, false, "<file>", "file: a mesh file", {}, {}, "file"},
}};

constexpr std::array<KeyInfo, 5> region_keys{{
    {"material", 0, false, "silicon", "the region's material (algaas takes its Al fraction)", materials, {}, {}},
    {"box", 0, true, "<lo hi> per dimension", "tensor meshes: the region's box", {}, {}, {}},
    {"polygon", 0, true, "<x y>, ...", "generate2d: the region's polygon", {}, {}, {}},
    {"seed", 0, true, "<x y z>", "generate3d with planes: a point inside the region", {}, {}, {}},
    {"band_to_band.direct_gap", 0, false, "off", "direct-gap band-to-band tunnelling", on_off, {}, {}},
}};

constexpr std::array<KeyInfo, 1> doping_keys{{
    {"profile", 0, true, "donor constant <peak>",
     "<donor|acceptor> <constant|gaussian|erfc> <peak> [axis=] [center=] [width=] [x=lo:hi] [lateral=]", dopants,
     shapes, {}},
}};

constexpr std::array<KeyInfo, 11> contact_keys{{
    {"kind", 0, false, "ohmic", "ohmic | gate | electrode | schottky", contact_kinds, {}, {}},
    {"patch", 0, false, "<boundary patch>", "the boundary patch the contact lies on", {}, {}, {}},
    {"box", 0, true, "<lo hi> per dimension", "only the patch's nodes inside the box", {}, {}, {}},
    {"oxide_thickness", 0, false, "<length>", "gate: the oxide's thickness", {}, {}, "gate"},
    {"oxide_permittivity", 0, false, "3.9", "gate: the oxide's relative permittivity", {}, {}, "gate"},
    {"fixed_charge", 0, false, "0", "gate: oxide fixed charge, cm^-2", {}, {}, "gate"},
    {"electrode", 0, false, "n_poly", "gate, electrode: n_poly | p_poly | metal", electrodes, {}, "gate electrode"},
    {"work_function", 0, false, "<eV>", "gate, electrode, schottky: a metal's work function, eV", {}, {}, "gate electrode schottky"},
    {"barrier", 0, false, "<eV>", "schottky: the electron barrier height, eV", {}, {}, "schottky"},
    {"richardson_n", 0, false, "0", "schottky: electron Richardson constant, A/(cm^2 K^2); 0 the material's", {}, {}, "schottky"},
    {"richardson_p", 0, false, "0", "schottky: hole Richardson constant, A/(cm^2 K^2); 0 the material's", {}, {}, "schottky"},
}};

constexpr std::array<KeyInfo, 8> interface_keys{{
    {"transport", 0, false, "drift_diffusion", "drift_diffusion | thermionic", transports, {}, {}},
    {"fixed_charge", 0, false, "0", "interface charge, cm^-2", {}, {}, {}},
    {"recombination_velocity_n", 0, false, "0", "surface recombination velocity, cm/s", {}, {}, {}},
    {"recombination_velocity_p", 0, false, "0", "surface recombination velocity, cm/s", {}, {}, {}},
    {"trap_level", 0, true, "acceptor <N_t> <E - E_i>", "a discrete trap level", dopants, {}, {}},
    {"trap_band", 0, true, "acceptor <D_it> <E_low> <E_high>", "a band of traps", dopants, {}, {}},
    {"thermal_velocity_n", 0, false, "<cm/s>", "traps: electron thermal velocity", {}, {}, {}},
    {"thermal_velocity_p", 0, false, "<cm/s>", "traps: hole thermal velocity", {}, {}, {}},
}};

constexpr std::array<KeyInfo, 4> thermal_contact_keys{{
    {"patch", 0, false, "<boundary patch>", "the boundary patch", {}, {}, {}},
    {"kind", 0, false, "isothermal", "isothermal | resistance", thermal_kinds, {}, {}},
    {"temperature", 0, false, "300", "K", {}, {}, {}},
    {"resistance", 0, false, "<K cm^2/W>", "resistance: the thermal resistance", {}, {}, {}},
}};

constexpr std::array<KeyInfo, 20> model_keys{{
    {"doping_mobility", 0, false, "on", "doping-dependent mobility", on_off, {}, {}},
    {"hole_mobility", 0, false, "automatic", "automatic | caughey_thomas | masetti", hole_mobility, {}, {}},
    {"electron_mobility", 0, false, "caughey_thomas", "caughey_thomas | masetti", electron_mobility, {}, {}},
    {"surface_mobility", 0, false, "off", "off | lombardi (inversion layers)", surface_mobility, {}, {}},
    {"quantum_confinement", 0, false, "off", "off | mlda (near insulators)", quantum_confinement, {}, {}},
    {"mlda_mass_n", 0, false, "0.98", "mlda: electron quantization mass, m0", {}, {}, {}},
    {"mlda_mass_p", 0, false, "0.49", "mlda: hole quantization mass, m0", {}, {}, {}},
    {"field_mobility", 0, false, "off", "high-field mobility (Canali)", on_off, {}, {}},
    {"srh", 0, false, "on", "Shockley-Read-Hall recombination", on_off, {}, {}},
    {"auger", 0, false, "on", "Auger recombination", on_off, {}, {}},
    {"radiative", 0, false, "on", "radiative recombination", on_off, {}, {}},
    {"bgn", 0, false, "on", "bandgap narrowing", on_off, {}, {}},
    {"fermi_dirac", 0, false, "off", "Fermi-Dirac statistics", on_off, {}, {}},
    {"incomplete_ionization", 0, false, "off", "incomplete dopant ionization", on_off, {}, {}},
    {"ionization_levels", 0, false, "fixed", "fixed | altermatt (with incomplete_ionization)", ionization_levels, {},
     {}},
    {"impact_ionization", 0, false, "off", "impact ionization", on_off, {}, {}},
    {"impact_current_resolution", 0, false, "0", "impact ionization: current resolution", {}, {}, {}},
    {"btbt_local", 0, false, "off", "local band-to-band tunnelling", on_off, {}, {}},
    {"btbt_nonlocal", 0, false, "off", "off | kane | direct_wkb", btbt_nonlocal, {}, {}},
    {"electrothermal", 0, false, "off", "lattice heat equation", on_off, {}, {}},
}};

constexpr std::array<KeyInfo, 26> solver_keys{{
    {"equations", 0, false, "drift_diffusion", "drift_diffusion | equilibrium_poisson", equations, {}, {}},
    {"max_iterations", 0, false, "100", "Newton iterations per point", {}, {}, {}},
    {"tol_update", 0, false, "1e-8", "convergence: largest scaled update", {}, {}, {}},
    {"max_update", 0, false, "5", "damping cap on the scaled update", {}, {}, {}},
    {"line_search", 0, false, "off", "off | residual | natural", line_search, {}, {}},
    {"line_search.max_backtracks", 0, false, "10", "line search: halvings", {}, {}, {}},
    {"line_search.sufficient_decrease", 0, false, "1e-4", "line search: decrease test", {}, {}, {}},
    {"max_cuts", 0, false, "0", "bias step cuts (0: off)", {}, {}, {}},
    {"growth", 0, false, "2", "sub-step growth after convergence", {}, {}, {}},
    {"cut_max_iterations", 0, false, "25", "Newton iterations per sub-step", {}, {}, {}},
    {"backend", 0, false, "eigen", "eigen | pardiso | gmres_ilut", backends, {}, {}},
    {"mkl", 0, false, "<path to mkl_rt>", "pardiso: the MKL runtime to load", {}, {}, {}},
    {"threads", 0, false, "1", "linear solver threads", {}, {}, {}},
    {"equilibrate", 0, false, "on", "scale the linear system", on_off, {}, {}},
    {"max_backward_error", 0, false, "1e-8", "linear solve acceptance", {}, {}, {}},
    {"max_refinement_steps", 0, false, "1", "iterative refinement steps", {}, {}, {}},
    {"min_pivot_ratio", 0, false, "1e-11", "smallest pivot ratio accepted", {}, {}, {}},
    {"restart", 0, false, "30", "gmres_ilut: restart length", {}, {}, {}},
    {"max_linear_iterations", 0, false, "1000", "gmres_ilut: iterations", {}, {}, {}},
    {"ilut_drop", 0, false, "1e-4", "gmres_ilut: drop tolerance", {}, {}, {}},
    {"ilut_fill", 0, false, "10", "gmres_ilut: fill factor", {}, {}, {}},
    {"linear_fallback", 0, false, "sparse_lu", "gmres_ilut: sparse_lu | none", fallbacks, {}, {}},
    {"Ns", 0, false, "<cm^-3>", "the density scale", {}, {}, {}},
    {"weighted_current", 0, false, "off", "weighted terminal current", on_off, {}, {}},
    {"rates", 0, false, "on", "store generation-recombination rates", on_off, {}, {}},
    {"minority_floor", 0, false, "off", "densities against eps of the majority", on_off, {}, {}},
}};

constexpr std::array<KeyInfo, 3> sweep_keys{{
    {"contact", 1, true, "<series>", "a contact's bias: a value, a list or first : last : step", {}, {}, {}},
    {"thermal", 1, true, "<series>", "a thermal contact's temperature, K", {}, {}, {}},
    {"restart", 0, false, "<checkpoint file>", "continue after the checkpoint's point", {}, {}, {}},
}};

constexpr std::array<KeyInfo, 12> transient_keys{{
    {"contact", 1, true, "constant(0)", "a contact's waveform", waveforms, {}, {}},
    {"t_end", 0, false, "<s>", "the run's end time (required)", {}, {}, {}},
    {"dt_initial", 0, false, "<s>", "the first time step", {}, {}, {}},
    {"dt_min", 0, false, "<s>", "the smallest time step", {}, {}, {}},
    {"dt_max", 0, false, "<s>", "the largest time step", {}, {}, {}},
    {"integrator", 0, false, "bdf2", "bdf2 | backward_euler", integrators, {}, {}},
    {"adaptive", 0, false, "on", "adaptive time steps", on_off, {}, {}},
    {"rtol", 0, false, "1e-3", "relative error tolerance", {}, {}, {}},
    {"density_ref", 0, false, "1e10", "densities below it are not error-controlled", {}, {}, {}},
    {"output_times", 0, false, "<series>", "times at which snapshots are kept", {}, {}, {}},
    {"fields_every_step", 0, false, "off", "a snapshot at every step", on_off, {}, {}},
    {"thermal", 1, true, "<K>", "a thermal contact's temperature, held", {}, {}, {}},
}};

constexpr std::array<KeyInfo, 4> small_signal_keys{{
    {"contact", 1, true, "<series>", "a contact's operating points", {}, {}, {}},
    {"thermal", 1, true, "<series>", "a thermal contact's temperature, K", {}, {}, {}},
    {"frequencies", 0, false, "<series>", "the frequencies, Hz (required)", {}, {}, {}},
    {"field_frequencies", 0, false, "<series>", "frequencies at which the fields are kept", {}, {}, {}},
}};

constexpr std::array<KeyInfo, 10> trace_keys{{
    {"contact", 0, false, "<name>", "the traced contact (required)", {}, {}, {}},
    {"end", 0, false, "<V>", "the traced contact's end bias (required)", {}, {}, {}},
    {"start", 1, true, "0", "a contact's bias at the start", {}, {}, {}},
    {"step", 0, false, "<V>", "the first arc-length step", {}, {}, {}},
    {"max_step", 0, false, "<V>", "the largest step", {}, {}, {}},
    {"min_step", 0, false, "<V>", "the smallest step", {}, {}, {}},
    {"current_limit", 0, false, "<A>", "stop when the current reaches it", {}, {}, {}},
    {"max_points", 0, false, "<n>", "the most points traced", {}, {}, {}},
    {"density_floor", 0, false, "<cm^-3>", "the density floor of the arc-length scale", {}, {}, {}},
    {"thermal", 1, true, "<K>", "a thermal contact's temperature, held", {}, {}, {}},
}};

constexpr std::array<KeyInfo, 7> output_keys{{
    {"sweep", 0, false, "<file>", "a sweep's or trace's results (.nitsweep)", {}, {}, {}},
    {"device", 0, false, "<file>", "the device (.nitdev)", {}, {}, {}},
    {"mesh", 0, false, "<file>", "the mesh", {}, {}, {}},
    {"checkpoint", 0, false, "<file>", "the last converged point", {}, {}, {}},
    {"table", 0, false, "<file>", "a CSV table of the terminal quantities", {}, {}, {}},
    {"transient", 0, false, "<file>", "a transient run's results (.nittran)", {}, {}, {}},
    {"small_signal", 0, false, "<file>", "a small-signal run's results (.nitac)", {}, {}, {}},
}};

constexpr std::array<SectionInfo, 15> sections{{
    {"device", 0, false, "the device's temperature", device_keys, {}},
    {"mesh", 0, false, "the mesh: tensor, generated or from a file", mesh_keys, "type"},
    {"region", 1, true, "a region: its material and where it is", region_keys, {}},
    {"doping", 0, false, "doping profiles", doping_keys, {}},
    {"contact", 1, true, "a contact on a boundary patch", contact_keys, "kind"},
    {"interface", 2, true, "the interface between two regions", interface_keys, {}},
    {"thermal_contact", 1, true, "a thermal contact", thermal_contact_keys, {}},
    {"models", 0, false, "physical models", model_keys, {}},
    {"solver", 0, false, "Newton and the linear solver", solver_keys, {}},
    {"equilibrium", 0, false, "the analysis: the equilibrium solve (no keys)", {}, {}},
    {"sweep", 0, false, "the analysis: a steady bias sweep", sweep_keys, {}},
    {"transient", 0, false, "the analysis: a transient run", transient_keys, {}},
    {"small_signal", 0, false, "the analysis: small-signal (AC) at operating points", small_signal_keys, {}},
    {"trace", 0, false, "the analysis: an arc-length trace", trace_keys, {}},
    {"output", 0, false, "files the run writes, relative to the deck", output_keys, {}},
}};

std::string_view trim(std::string_view s) {
    const auto b = s.find_first_not_of(" \t\r");
    if (b == std::string_view::npos) return {};
    const auto e = s.find_last_not_of(" \t\r");
    return s.substr(b, e - b + 1);
}

std::string_view first_word(std::string_view s) {
    s = trim(s);
    return s.substr(0, s.find_first_of(" \t"));
}

std::vector<std::string_view> words_of(std::string_view s) {
    std::vector<std::string_view> out;
    std::size_t i = 0;
    while (i < s.size()) {
        while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == ',')) ++i;
        const std::size_t b = i;
        while (i < s.size() && s[i] != ' ' && s[i] != '\t' && s[i] != ',') ++i;
        if (i > b) out.push_back(s.substr(b, i - b));
    }
    return out;
}

bool has_context(const KeyInfo& k, std::string_view context) {
    return std::ranges::contains(words_of(k.contexts), context);
}

std::vector<std::string_view> lines_of(std::string_view text) {
    std::vector<std::string_view> out;
    std::size_t start = 0;
    while (start <= text.size()) {
        const auto end = text.find('\n', start);
        auto l = text.substr(start, end == std::string_view::npos ? std::string_view::npos : end - start);
        if (!l.empty() && l.back() == '\r') l.remove_suffix(1);
        out.push_back(l);
        if (end == std::string_view::npos) break;
        start = end + 1;
    }
    return out;
}

// The part of a line before its comment ('#' outside double quotes).
std::string_view uncommented(std::string_view line) {
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        if (line[i] == '"') quoted = !quoted;
        else if (line[i] == '#' && !quoted) return line.substr(0, i);
    }
    return line;
}

bool continues(std::string_view line) {
    const auto t = trim(uncommented(line));
    return !t.empty() && t.back() == '\\';
}

bool identifier_char(char c) {
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_' || c == '.' ||
           c == '-';
}

// The section a line belongs to (the last header at or above it), with its context value.
struct Where {
    const SectionLocation* section = nullptr;
    std::string context;
};

Where where(const std::vector<SectionLocation>& found, int line) {
    Where w;
    for (const SectionLocation& s : found) {
        if (s.line <= line) w.section = &s;
    }
    if (!w.section) return w;
    if (const SectionInfo* info = find_section_info(w.section->kind); info && !info->context_key.empty()) {
        for (const auto& e : w.section->entries) {
            if (e.key == info->context_key) w.context = std::string(first_word(e.value));
        }
    }
    return w;
}

std::vector<std::string_view> key_names(const SectionInfo& s, std::string_view context) {
    std::vector<std::string_view> out;
    for (const KeyInfo* k : section_keys(s, context)) {
        if (std::ranges::find(out, k->key) == out.end()) out.push_back(k->key);
    }
    return out;
}

std::string did_you_mean(std::string_view name) { return "did you mean '" + std::string(name) + "'?"; }

bool starts_with(std::string_view s, std::string_view prefix) { return s.substr(0, prefix.size()) == prefix; }

}  // namespace

std::span<const SectionInfo> deck_vocabulary() noexcept { return sections; }

const SectionInfo* find_section_info(std::string_view kind) noexcept {
    for (const SectionInfo& s : sections) {
        if (s.kind == kind) return &s;
    }
    return nullptr;
}

std::vector<const KeyInfo*> section_keys(const SectionInfo& s, std::string_view context) {
    // A context that is not one of the selecting key's choices (empty, or a misspelled type) selects
    // every key.
    const auto selector = std::ranges::find(s.keys, s.context_key, &KeyInfo::key);
    const bool known = !s.context_key.empty() && selector != s.keys.end() &&
                       std::ranges::contains(selector->choices, context);
    std::vector<const KeyInfo*> out;
    for (const KeyInfo& k : s.keys) {
        if (!known || k.contexts.empty() || has_context(k, context)) out.push_back(&k);
    }
    return out;
}

const KeyInfo* find_key_info(const SectionInfo& s, std::string_view key, std::string_view context) {
    const auto words = words_of(key);
    if (words.empty()) return nullptr;
    const KeyInfo* any = nullptr;
    for (const KeyInfo* k : section_keys(s, context)) {
        if (k->key != words.front()) continue;
        if (static_cast<std::size_t>(k->names) + 1 == words.size()) return k;
        if (!any) any = k;
    }
    return any;
}

std::size_t edit_distance(std::string_view a, std::string_view b) {
    std::vector<std::size_t> row(b.size() + 1);
    std::iota(row.begin(), row.end(), std::size_t{0});
    for (std::size_t i = 1; i <= a.size(); ++i) {
        std::size_t diagonal = row[0];
        row[0] = i;
        for (std::size_t j = 1; j <= b.size(); ++j) {
            const std::size_t up = row[j];
            row[j] = std::min({row[j] + 1, row[j - 1] + 1, diagonal + (a[i - 1] == b[j - 1] ? 0 : 1)});
            diagonal = up;
        }
    }
    return row[b.size()];
}

std::optional<std::string> nearest(std::string_view word, std::span<const std::string_view> names) {
    if (word.empty() || std::ranges::find(names, word) != names.end()) return std::nullopt;
    const std::size_t limit = std::min<std::size_t>(2, std::max<std::size_t>(1, word.size() / 3));
    std::optional<std::string_view> best;
    std::size_t best_distance = limit + 1;
    for (const std::string_view n : names) {
        if (n.empty()) continue;
        const std::size_t d = edit_distance(word, n);
        if (d < best_distance) {
            best = n;
            best_distance = d;
        }
    }
    if (!best) return std::nullopt;
    return std::string(*best);
}

std::optional<std::string> suggestion(std::string_view text, const Diagnostic& d) {
    const auto lines = lines_of(text);
    if (d.line < 1 || static_cast<std::size_t>(d.line) > lines.size()) return std::nullopt;
    const std::string_view line = trim(uncommented(lines[static_cast<std::size_t>(d.line) - 1]));
    if (line.empty()) return std::nullopt;
    if (line.front() == '[') {
        const auto close = line.find(']');
        const auto kind = first_word(line.substr(1, close == std::string_view::npos ? std::string_view::npos : close - 1));
        if (find_section_info(kind)) return std::nullopt;
        std::vector<std::string_view> kinds;
        for (const SectionInfo& s : sections) kinds.push_back(s.kind);
        if (auto n = nearest(kind, kinds)) return did_you_mean(*n);
        return std::nullopt;
    }
    const auto eq = line.find('=');
    if (eq == std::string_view::npos) return std::nullopt;
    const auto found = find_sections(text);
    const Where w = where(found, d.line);
    if (!w.section) return std::nullopt;
    const SectionInfo* info = find_section_info(w.section->kind);
    if (!info) return std::nullopt;
    const std::string_view key = trim(line.substr(0, eq));
    const KeyInfo* k = find_key_info(*info, key, w.context);
    if (!k) {
        if (auto n = nearest(first_word(key), key_names(*info, w.context))) return did_you_mean(*n);
        return std::nullopt;
    }
    const auto value = words_of(trim(line.substr(eq + 1)));
    // A waveform is a call: its name is the text before '('.
    const auto name_of = [](std::string_view word) { return word.substr(0, word.find('(')); };
    if (!k->choices.empty() && !value.empty()) {
        if (auto n = nearest(name_of(value[0]), k->choices)) return did_you_mean(*n);
    }
    if (!k->choices2.empty() && value.size() > 1) {
        if (auto n = nearest(value[1], k->choices2)) return did_you_mean(*n);
    }
    return std::nullopt;
}

std::optional<Completion> complete(std::string_view text, int line, std::size_t column) {
    const auto lines = lines_of(text);
    if (line < 1 || static_cast<std::size_t>(line) > lines.size()) return std::nullopt;
    const std::string_view l = lines[static_cast<std::size_t>(line) - 1];
    column = std::min(column, l.size());
    if (line > 1 && continues(lines[static_cast<std::size_t>(line) - 2])) return std::nullopt;
    const std::string_view before = l.substr(0, column);
    if (uncommented(before).size() < before.size()) return std::nullopt;  // in a comment
    // The word at the caret: identifier characters on both sides.
    std::size_t from = column, to = column;
    while (from > 0 && identifier_char(l[from - 1])) --from;
    while (to < l.size() && identifier_char(l[to])) ++to;
    const std::string_view typed = l.substr(from, column - from);
    Completion c{line, from, to, {}};

    const std::string_view lead = trim(before);
    if (!lead.empty() && lead.front() == '[') {
        // A header: only its kind, the first word after '['.
        if (lead.find(']') != std::string_view::npos) return std::nullopt;
        if (trim(lead.substr(1)) != typed) return std::nullopt;
        const bool closed = l.find(']', to) != std::string_view::npos;
        for (const SectionInfo& s : sections) {
            if (!starts_with(s.kind, typed)) continue;
            std::string insert(s.kind);
            if (s.names > 0) insert += " ";
            else if (!closed) insert += "]";
            c.items.push_back({std::string(s.kind), std::move(insert), std::string(s.help)});
        }
        return c;
    }
    const auto found = find_sections(text);
    const Where w = where(found, line);
    if (!w.section || w.section->line == line) return std::nullopt;
    const SectionInfo* info = find_section_info(w.section->kind);
    if (!info) return std::nullopt;

    const auto eq = before.find('=');
    if (eq == std::string_view::npos) {
        // An entry's key: the first word of the line.
        if (lead != typed) return std::nullopt;
        for (const KeyInfo* k : section_keys(*info, w.context)) {
            if (!starts_with(k->key, typed)) continue;
            const bool present = k->names == 0 && !k->repeatable &&
                                 std::ranges::any_of(w.section->entries, [&](const SectionLocation::Entry& e) {
                                     return e.line != line && e.key == k->key;
                                 });
            if (present) continue;
            const bool listed = std::ranges::any_of(c.items, [&](const Candidate& o) { return o.label == k->key; });
            if (listed) continue;
            std::string insert(k->key);
            insert += k->names > 0 ? " " : " = ";
            c.items.push_back({std::string(k->key), std::move(insert), std::string(k->help)});
        }
        return c;
    }
    // A value: its first word's choices, or the second word's.
    const KeyInfo* k = find_key_info(*info, trim(before.substr(0, eq)), w.context);
    if (!k) return std::nullopt;
    const auto words = words_of(before.substr(eq + 1));
    const bool in_word = !typed.empty() || (column > 0 && identifier_char(l[column - 1]));
    const std::size_t index = words.empty() ? 0 : words.size() - (in_word ? 1 : 0);
    const Words choices = index == 0 ? k->choices : index == 1 ? k->choices2 : Words{};
    for (const std::string_view choice : choices) {
        if (!starts_with(choice, typed)) continue;
        std::string insert(choice);
        if (k->choices.data() == waveforms.data()) insert += "(";
        c.items.push_back({std::string(choice), std::move(insert), std::string(k->help)});
    }
    return c;
}

}  // namespace NiTCAD::app
