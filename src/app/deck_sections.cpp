// The menus' deck sections (deck_sections.hpp). The defaults are those deck/build.hpp documents; the
// WG7 test checks that inserting each one leaves the plan unchanged.
#include "NiTCAD/app/deck_sections.hpp"

#include <algorithm>
#include <array>

namespace NiTCAD::app {

namespace {

constexpr std::array<KeySpec, 1> device_keys{{{"temperature", "300", "lattice temperature, K"}}};

// [mesh]: only a generated mesh's options with defaults; the size keys have none.
constexpr std::array<KeySpec, 3> mesh_keys{{{"min_angle", "20.7", "generate2d: smallest triangle angle, degrees"},
                                            {"radius_edge", "2", "generate3d: largest radius-edge ratio"},
                                            {"max_nodes", "2000000", "generated meshes: node limit"}}};

constexpr std::array<KeySpec, 17> model_keys{{
    {"doping_mobility", "on", "doping-dependent mobility"},
    {"hole_mobility", "automatic", "automatic | caughey_thomas | masetti"},
    {"electron_mobility", "caughey_thomas", "caughey_thomas | masetti"},
    {"surface_mobility", "off", "off | lombardi (inversion layers)"},
    {"field_mobility", "off", "high-field mobility (Canali)"},
    {"srh", "on", "Shockley-Read-Hall recombination"},
    {"auger", "on", "Auger recombination"},
    {"radiative", "on", "radiative recombination"},
    {"bgn", "on", "bandgap narrowing"},
    {"fermi_dirac", "off", "Fermi-Dirac statistics"},
    {"incomplete_ionization", "off", "incomplete dopant ionization"},
    {"ionization_levels", "fixed", "fixed | altermatt (with incomplete_ionization)"},
    {"impact_ionization", "off", "impact ionization"},
    {"impact_current_resolution", "0", "impact ionization: current resolution"},
    {"btbt_local", "off", "local band-to-band tunnelling"},
    {"btbt_nonlocal", "off", "off | kane | direct_wkb"},
    {"electrothermal", "off", "lattice heat equation"},
}};

constexpr std::array<KeySpec, 23> solver_keys{{
    {"equations", "drift_diffusion", "drift_diffusion | equilibrium_poisson"},
    {"max_iterations", "100", "Newton iterations per point"},
    {"tol_update", "1e-8", "convergence: largest scaled update"},
    {"max_update", "5", "damping cap on the scaled update"},
    {"line_search", "off", "off | residual | natural"},
    {"line_search.max_backtracks", "10", "line search: halvings"},
    {"line_search.sufficient_decrease", "1e-4", "line search: decrease test"},
    {"max_cuts", "0", "bias step cuts (0: off)"},
    {"growth", "2", "sub-step growth after convergence"},
    {"cut_max_iterations", "25", "Newton iterations per sub-step"},
    {"backend", "eigen", "eigen | pardiso | gmres_ilut"},
    {"threads", "1", "linear solver threads"},
    {"equilibrate", "on", "scale the linear system"},
    {"max_backward_error", "1e-8", "linear solve acceptance"},
    {"max_refinement_steps", "1", "iterative refinement steps"},
    {"min_pivot_ratio", "1e-11", "smallest pivot ratio accepted"},
    {"restart", "30", "gmres_ilut: restart length"},
    {"max_linear_iterations", "1000", "gmres_ilut: iterations"},
    {"ilut_drop", "1e-4", "gmres_ilut: drop tolerance"},
    {"ilut_fill", "10", "gmres_ilut: fill factor"},
    {"linear_fallback", "sparse_lu", "gmres_ilut: sparse_lu | none"},
    {"weighted_current", "off", "weighted terminal current"},
    {"rates", "on", "store generation-recombination rates"},
}};

constexpr std::array<KeySpec, 5> transient_keys{{{"integrator", "bdf2", "bdf2 | backward_euler"},
                                                 {"adaptive", "on", "adaptive time steps"},
                                                 {"rtol", "1e-3", "relative error tolerance"},
                                                 {"density_ref", "1e10", "densities below it are not error-controlled"},
                                                 {"fields_every_step", "off", "a snapshot at every step"}}};

constexpr std::array<SectionSpec, 6> device_sections{{{"device", "[device]", device_keys},
                                                      {"region", "[region]", {}},
                                                      {"doping", "[doping]", {}},
                                                      {"contact", "[contact]", {}},
                                                      {"interface", "[interface]", {}},
                                                      {"thermal_contact", "[thermal_contact]", {}}}};
constexpr std::array<SectionSpec, 1> mesh_sections{{{"mesh", "[mesh]", mesh_keys}}};
constexpr std::array<SectionSpec, 1> physics_sections{{{"models", "[models]", model_keys}}};
constexpr std::array<SectionSpec, 8> simulation_sections{{{"solver", "[solver]", solver_keys},
                                                          {"equilibrium", "[equilibrium]", {}},
                                                          {"sweep", "[sweep]", {}},
                                                          {"transient", "[transient]", transient_keys},
                                                          {"small_signal", "[small_signal]", {}},
                                                          {"trace", "[trace]", {}},
                                                          {"output", "[output]", {}},
                                                          {"thermal_contact", "[thermal_contact]", {}}}};

bool is_analysis(std::string_view kind) {
    return kind == "equilibrium" || kind == "sweep" || kind == "transient" || kind == "small_signal" || kind == "trace";
}

std::string_view trim(std::string_view s) {
    const auto b = s.find_first_not_of(" \t\r");
    if (b == std::string_view::npos) return {};
    const auto e = s.find_last_not_of(" \t\r");
    return s.substr(b, e - b + 1);
}

// Runs of blanks collapsed to one, as deck/syntax.hpp collapses names and keys.
std::string collapse(std::string_view s) {
    std::string out;
    bool blank = false;
    for (const char c : trim(s)) {
        if (c == ' ' || c == '\t') {
            blank = true;
            continue;
        }
        if (blank && !out.empty()) out += ' ';
        blank = false;
        out += c;
    }
    return out;
}

std::vector<std::string_view> lines_of(std::string_view text) {
    std::vector<std::string_view> out;
    std::size_t start = 0;
    while (start <= text.size()) {
        const auto end = text.find('\n', start);
        out.push_back(text.substr(start, end == std::string_view::npos ? std::string_view::npos : end - start));
        if (end == std::string_view::npos) break;
        start = end + 1;
    }
    if (!out.empty() && out.back().empty() && !text.empty()) out.pop_back();  // a final newline
    return out;
}

}  // namespace

std::string_view to_string(DeckMenu m) noexcept {
    switch (m) {
    case DeckMenu::device: return "Device";
    case DeckMenu::mesh: return "Mesh";
    case DeckMenu::physics: return "Physics";
    case DeckMenu::simulation: return "Simulation";
    }
    return "?";
}

std::span<const SectionSpec> menu_sections(DeckMenu m) noexcept {
    switch (m) {
    case DeckMenu::device: return device_sections;
    case DeckMenu::mesh: return mesh_sections;
    case DeckMenu::physics: return physics_sections;
    case DeckMenu::simulation: return simulation_sections;
    }
    return {};
}

const SectionSpec* find_section_spec(std::string_view kind) noexcept {
    for (const DeckMenu m : {DeckMenu::device, DeckMenu::mesh, DeckMenu::physics, DeckMenu::simulation}) {
        for (const SectionSpec& s : menu_sections(m)) {
            if (s.kind == kind) return &s;
        }
    }
    return nullptr;
}

std::vector<SectionLocation> find_sections(std::string_view text) {
    std::vector<SectionLocation> out;
    const auto lines = lines_of(text);
    bool continued = false;  // the previous line ended in a backslash: this one is part of its value
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const int number = static_cast<int>(i) + 1;
        std::string_view line = lines[i];
        // The comment, outside double quotes.
        bool quoted = false;
        for (std::size_t k = 0; k < line.size(); ++k) {
            if (line[k] == '"') quoted = !quoted;
            if (line[k] == '#' && !quoted) {
                line = line.substr(0, k);
                break;
            }
        }
        const std::string_view t = trim(line);
        const bool was_continued = continued;
        continued = !t.empty() && t.back() == '\\';
        if (was_continued) {
            if (!out.empty()) out.back().last_line = number;
            continue;
        }
        if (t.empty()) continue;
        if (t.front() == '[') {
            const auto close = t.find(']');
            const std::string inside = collapse(t.substr(1, close == std::string_view::npos ? std::string_view::npos : close - 1));
            SectionLocation s;
            s.line = s.last_line = number;
            const auto space = inside.find(' ');
            s.kind = inside.substr(0, space);
            if (space != std::string::npos) s.name = inside.substr(space + 1);
            out.push_back(std::move(s));
            continued = false;
            continue;
        }
        if (out.empty()) continue;
        const auto eq = t.find('=');
        if (eq != std::string_view::npos) {
            out.back().entries.push_back({collapse(t.substr(0, eq)), std::string(trim(t.substr(eq + 1))), number});
        }
        out.back().last_line = number;
    }
    return out;
}

std::optional<SectionLocation> find_section(std::string_view text, std::string_view kind) {
    for (auto& s : find_sections(text)) {
        if (s.kind == kind) return std::move(s);
    }
    return std::nullopt;
}

std::vector<KeySpec> insertable_keys(const SectionLocation& section) {
    std::vector<KeySpec> out;
    const SectionSpec* spec = find_section_spec(section.kind);
    if (!spec) return out;
    // A mesh's options depend on its type.
    std::string type;
    for (const auto& e : section.entries) {
        if (e.key == "type") type = e.value;
    }
    for (const KeySpec& k : spec->keys) {
        if (std::ranges::any_of(section.entries, [&](const auto& e) { return e.key == k.key; })) continue;
        if (section.kind == "mesh") {
            if (type != "generate2d" && type != "generate3d") continue;
            if (k.key == "min_angle" && type != "generate2d") continue;
            if (k.key == "radius_edge" && type != "generate3d") continue;
        }
        out.push_back(k);
    }
    return out;
}

std::optional<Insertion> insert_key(std::string_view text, std::string_view kind, const KeySpec& key) {
    const auto section = find_section(text, kind);
    const std::string entry = std::string(key.key) + " = " + std::string(key.value);
    if (!section) {
        if (is_analysis(kind)) return std::nullopt;
        Insertion out{std::string(text), 0};
        if (!out.text.empty() && out.text.back() != '\n') out.text += '\n';
        out.text += "\n[" + std::string(kind) + "]\n" + entry + "\n";
        out.line = static_cast<int>(std::ranges::count(out.text, '\n'));
        return out;
    }
    // After the section's last entry line: find where that line ends.
    std::size_t at = 0;
    for (int line = 1; line <= section->last_line; ++line) {
        const auto nl = text.find('\n', at);
        if (nl == std::string_view::npos) {
            at = text.size();
            break;
        }
        at = nl + 1;
    }
    Insertion out;
    out.text = std::string(text.substr(0, at));
    if (!out.text.empty() && out.text.back() != '\n') out.text += '\n';
    out.text += entry + "\n";
    out.text += text.substr(at);
    out.line = section->last_line + 1;
    return out;
}

}  // namespace NiTCAD::app
