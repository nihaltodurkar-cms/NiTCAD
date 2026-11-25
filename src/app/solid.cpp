// The 3D solid modeller's model (solid.hpp; DECISIONS.md Unit 26l, SD3).
#include "NiTCAD/app/solid.hpp"

#include "NiTCAD/app/dimensions.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <numbers>
#include <ranges>
#include <set>

#include "deck_text.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/mesh/generate.hpp"

namespace NiTCAD::app::solid {

namespace {

using namespace deck_text;
using geo::format_length;
using geo::parse_length;

using Params = std::vector<deck::Parameter>;

// One coordinate: a value, or a reference resolved by `params` (Unit 26m MT8); its reference ("" for a value) appended
// to `refs`.
std::optional<Coord> coord(std::string_view w, const Params& params, std::vector<std::string>& refs) {
    if (w.find('$') == std::string_view::npos) {
        refs.emplace_back();
        return parse_length(w);
    }
    const auto t = deck::substitute(w, params);
    if (!t) return std::nullopt;
    refs.emplace_back(w);
    return parse_length(*t);
}

// No reference among them: none kept (Solid::refs empty).
void tidy(std::vector<std::string>& refs) {
    if (std::ranges::all_of(refs, [](const std::string& r) { return r.empty(); })) refs.clear();
}

std::optional<std::vector<Coord>> lengths(std::string_view value, std::size_t n, const Params& params, std::vector<std::string>& refs) {
    const auto w = deck::value_words(value);
    if (w.size() != n) return std::nullopt;
    std::vector<Coord> out;
    for (const auto t : w) {
        const auto v = coord(t, params, refs);
        if (!v) return std::nullopt;
        out.push_back(*v);
    }
    return out;
}

// A region box or the domain: x0 x1 y0 y1 z0 z1.
std::optional<Box3> parse_box3(std::string_view value, const Params& params, std::vector<std::string>& refs) {
    const auto c = lengths(value, 6, params, refs);
    if (!c || (*c)[0] > (*c)[1] || (*c)[2] > (*c)[3] || (*c)[4] > (*c)[5]) return std::nullopt;
    tidy(refs);
    return Box3{{(*c)[0], (*c)[2], (*c)[4]}, {(*c)[1], (*c)[3], (*c)[5]}};
}

std::optional<Box3> parse_box3(std::string_view value, const Params& params) {
    std::vector<std::string> refs;
    return parse_box3(value, params, refs);
}

// A patch box: x0 y0 z0 x1 y1 z1 (its references read, not kept: a patch is written as values).
std::optional<Box3> parse_patch(std::string_view value, const Params& params) {
    std::vector<std::string> refs;
    const auto c = lengths(value, 6, params, refs);
    if (!c || (*c)[0] > (*c)[3] || (*c)[1] > (*c)[4] || (*c)[2] > (*c)[5]) return std::nullopt;
    return Box3{{(*c)[0], (*c)[1], (*c)[2]}, {(*c)[3], (*c)[4], (*c)[5]}};
}

std::optional<std::vector<Pt>> parse_points(std::string_view value, const Params& params, std::vector<std::string>& refs) {
    const auto w = deck::value_words(value);
    if (w.size() < 6 || w.size() % 2 != 0) return std::nullopt;
    std::vector<Pt> out;
    for (std::size_t k = 0; k < w.size(); k += 2) {
        const auto u = coord(w[k], params, refs);
        const auto v = coord(w[k + 1], params, refs);
        if (!u || !v) return std::nullopt;
        out.push_back({*u, *v});
    }
    return out;
}

std::optional<std::vector<Pt>> parse_points(std::string_view value) {
    std::vector<std::string> refs;
    return parse_points(value, {}, refs);
}

std::optional<int> axis_of(std::string_view a) {
    if (a == "x") return 0;
    if (a == "y") return 1;
    if (a == "z") return 2;
    return std::nullopt;
}

std::optional<Solid> parse_solid(const deck::Entry& e, const Params& params) {
    if (e.key == "box") {
        std::vector<std::string> refs;
        const auto b = parse_box3(e.value, params, refs);
        if (!b) return std::nullopt;
        Solid s{SolidKind::box, *b, 2, 0, 0, {}, {}};
        s.refs = std::move(refs);
        return s;
    }
    if (e.key == "prism") {
        // <axis> <from> <to>, u v, u v, ...
        const std::string_view text = e.value;
        const std::size_t comma = text.find(',');
        if (comma == std::string_view::npos) return std::nullopt;
        const auto head = deck::value_words(text.substr(0, comma));
        if (head.size() != 3) return std::nullopt;
        const auto axis = axis_of(head[0]);
        Solid s;
        const auto from = coord(head[1], params, s.refs);
        const auto to = coord(head[2], params, s.refs);
        auto p = parse_points(text.substr(comma + 1), params, s.refs);
        if (!axis || !from || !to || !p) return std::nullopt;
        s.kind = SolidKind::prism;
        s.axis = *axis;
        s.from = *from;
        s.to = *to;
        s.polygon = std::move(*p);
        tidy(s.refs);
        return s;
    }
    return std::nullopt;
}

bool is_solid_key(std::string_view k) { return k == "box" || k == "prism"; }

constexpr const char* axis_name[] = {"x", "y", "z"};

// The model's parameters as the deck layer resolves them (for deck::substitute).
Params deck_params(const Model& m) {
    Params out;
    for (const Param& p : m.parameters) out.push_back({p.name, p.text, p.length ? static_cast<double>(p.value) * geo::quantum_cm : p.number, p.length, 0});
    return out;
}

// What a reference gives, in quanta; none if it does not resolve or is not a length.
std::optional<Coord> resolved(std::string_view ref, const Params& params) {
    const auto t = deck::substitute(ref, params);
    if (!t) return std::nullopt;
    return parse_length(*t);
}

// A coordinate as written: its reference when it still gives the value, else the value.
std::string coord_text(Coord v, const std::vector<std::string>& refs, std::size_t k, const Params& params) {
    if (k < refs.size() && !refs[k].empty() && resolved(refs[k], params) == v) return refs[k];
    return format_length(v);
}

std::string box_text(const Box3& b, const std::vector<std::string>& refs = {}, const Params& params = {}) {  // x0 x1 y0 y1 z0 z1
    const Coord c[6] = {b.lo.x, b.hi.x, b.lo.y, b.hi.y, b.lo.z, b.hi.z};
    std::string out;
    for (std::size_t k = 0; k < 6; ++k) out += (k ? " " : "") + coord_text(c[k], refs, k, params);
    return out;
}

std::string patch_text(const Box3& b) {  // x0 y0 z0  x1 y1 z1
    return format_length(b.lo.x) + " " + format_length(b.lo.y) + " " + format_length(b.lo.z) + "  " + format_length(b.hi.x) + " " +
           format_length(b.hi.y) + " " + format_length(b.hi.z);
}

std::string solid_line(const Solid& s, const Params& params = {}) {
    if (s.kind == SolidKind::box) return "box = " + box_text(s.box, s.refs, params);
    std::string out = std::string("prism = ") + axis_name[s.axis] + " " + coord_text(s.from, s.refs, 0, params) + " " +
                      coord_text(s.to, s.refs, 1, params);
    for (std::size_t j = 0; j < s.polygon.size(); ++j) {
        out += ", " + coord_text(s.polygon[j].x, s.refs, 2 + 2 * j, params) + " " + coord_text(s.polygon[j].y, s.refs, 3 + 2 * j, params);
    }
    return out;
}

// Twice the polygon's signed area (exact for coordinates within the limit and a handful of vertices).
std::int64_t twice_area(const std::vector<Pt>& p) {
    std::int64_t a = 0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        const Pt q = p[i], r = p[(i + 1) % p.size()];
        a += q.x * r.y - r.x * q.y;
    }
    return a;
}

bool identifier(std::string_view s) {
    if (s.empty() || !std::isalpha(static_cast<unsigned char>(s.front()))) return false;
    return std::ranges::all_of(s, [](char c) { return std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '.' || c == '-'; });
}

bool within(Coord c) { return c >= -coord_limit && c <= coord_limit; }
bool within(P3 p) { return within(p.x) && within(p.y) && within(p.z); }

bool box_ok(const Box3& b) { return b.lo.x < b.hi.x && b.lo.y < b.hi.y && b.lo.z < b.hi.z && within(b.lo) && within(b.hi); }

}  // namespace

Box3 normalized(P3 a, P3 b) noexcept {
    return {{std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)}, {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)}};
}

P3 from_uv(int axis, Pt uv, Coord t) noexcept {
    P3 p;
    p[axis] = t;
    p[(axis + 1) % 3] = uv.x;
    p[(axis + 2) % 3] = uv.y;
    return p;
}

Pt uv_of(int axis, P3 p) noexcept { return {p[(axis + 1) % 3], p[(axis + 2) % 3]}; }

Box3 Solid::bounds() const noexcept {
    if (kind == SolidKind::box) return box;
    Pt lo = polygon.front(), hi = polygon.front();
    for (const Pt& p : polygon) {
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y)};
        hi = {std::max(hi.x, p.x), std::max(hi.y, p.y)};
    }
    return normalized(from_uv(axis, lo, from), from_uv(axis, hi, to));
}

std::vector<P3> Solid::vertices() const {
    std::vector<P3> out;
    if (kind == SolidKind::box) {
        for (int k = 0; k < 8; ++k) out.push_back({(k & 1) ? box.hi.x : box.lo.x, (k & 2) ? box.hi.y : box.lo.y, (k & 4) ? box.hi.z : box.lo.z});
        return out;
    }
    for (const Coord t : {from, to}) {
        for (const Pt& p : polygon) out.push_back(from_uv(axis, p, t));
    }
    return out;
}

std::vector<std::pair<int, int>> Solid::edges() const {
    std::vector<std::pair<int, int>> out;
    if (kind == SolidKind::box) {
        for (int k = 0; k < 8; ++k) {
            for (const int bit : {1, 2, 4}) {
                if (!(k & bit)) out.emplace_back(k, k | bit);
            }
        }
        return out;
    }
    const int n = static_cast<int>(polygon.size());
    for (int k = 0; k < n; ++k) {
        out.emplace_back(k, (k + 1) % n);          // the cap at `from`
        out.emplace_back(n + k, n + (k + 1) % n);  // at `to`
        out.emplace_back(k, n + k);                // along the axis
    }
    return out;
}

int Model::background() const noexcept {
    for (std::size_t k = 0; k < regions.size(); ++k) {
        if (regions[k].solids.empty()) return static_cast<int>(k);
    }
    return -1;
}

// --- reading -------------------------------------------------------------------------------------------------------

std::expected<Model, std::string> read_model(std::string_view text) {
    const auto d = deck::parse(text);
    if (!d) return std::unexpected(d.error().message);
    Model m;
    const deck::Section* mesh = nullptr;
    for (const auto& s : d->sections) {
        if (s.kind == "mesh") mesh = &s;
    }
    const auto refuse = [&](std::string why) {
        m.mesh = MeshKind::unsupported;
        m.why = std::move(why);
        m.regions.clear();
        m.contacts.clear();
        return m;
    };
    const auto at_line = [](const deck::Entry& e, std::string what) { return "line " + std::to_string(e.line) + ": " + what; };
    if (!mesh) return refuse("the deck has no [mesh] section");
    // Unit 26m MT8: the parameters, each kept with its expression as written.
    const auto resolved_params = deck::parameters(*d);
    if (!resolved_params) return refuse(resolved_params.error().message);
    const Params& params = *resolved_params;
    for (const auto& s : d->sections) {
        if (s.kind != "parameters") continue;
        for (std::size_t k = 0; k < s.entries.size(); ++k) {
            const auto it = std::ranges::find(params, s.entries[k].key, &deck::Parameter::name);
            Param p{s.entries[k].key, std::string(trim(s.entries[k].value)), it->text, it->length, 0, it->length ? 0.0 : it->value};
            if (it->length) {
                const auto q = parse_length(it->text);
                if (!q) return refuse(at_line(s.entries[k], "parameter '" + p.name + "' is not a length the modeller can draw"));
                p.value = *q;
            }
            m.parameters.push_back(std::move(p));
        }
    }
    const deck::Entry* type = find_entry(*mesh, "type");
    const std::string kind = type ? std::string(trim(type->value)) : std::string("tensor");
    if (kind == "generate3d") {
        const deck::Entry* l = find_entry(*mesh, "layered");
        if (!l || !on_off(l->value).value_or(false)) {
            return refuse("a generate3d deck without layers (a cuboid cut by planes) is not drawn: the modeller draws layered ones");
        }
        m.mesh = MeshKind::layered;
    } else if (kind == "extrude3d") {
        m.mesh = MeshKind::extruded;
    } else if (kind == "tensor") {
        if (find_entry(*mesh, "z")) return refuse("a 3D tensor deck is not drawn in the modeller (layered generate3d and extrude3d decks are)");
        return refuse("a 1D or 2D device is drawn in the 2D editor");
    } else if (kind == "generate2d") {
        return refuse("a 2D device is drawn in the 2D editor");
    } else {
        return refuse("a " + kind + " mesh is not drawn in the modeller (layered generate3d and extrude3d decks are)");
    }

    std::map<std::string, Box3, std::less<>> patches;
    int axis = 0;
    std::vector<Coord> stations;
    std::map<std::string, std::vector<Pt>, std::less<>> cells;
    std::vector<std::string> cell_order;
    bool domain = false;
    for (const auto& e : mesh->entries) {
        if (e.key.starts_with("patch ")) {
            const auto w = words(e.key);
            const auto b = parse_patch(e.value, params);
            if (w.size() != 2 || !b) return refuse(at_line(e, "a patch box the modeller cannot read"));
            patches[std::string(w[1])] = *b;
        } else if (m.mesh == MeshKind::layered && e.key == "domain") {
            const auto b = parse_box3(e.value, params, m.domain_refs);
            if (!b) return refuse(at_line(e, "a domain the modeller cannot read"));
            m.domain = *b;
            domain = true;
        } else if (m.mesh == MeshKind::extruded && e.key == "axis") {
            const auto a = axis_of(trim(e.value));
            if (!a) return refuse(at_line(e, "an axis the modeller cannot read"));
            axis = *a;
        } else if (m.mesh == MeshKind::extruded && e.key == "stations") {
            for (const auto t : words(e.value)) {
                const auto v = parse_length(t);
                if (!v) return refuse(at_line(e, "a station the modeller cannot read"));
                stations.push_back(*v);
            }
        } else if (m.mesh == MeshKind::extruded && e.key.starts_with("cell ")) {
            const auto w = words(e.key);
            auto p = parse_points(e.value);
            if (w.size() != 2 || !p) return refuse(at_line(e, "a cell the modeller cannot read"));
            cells[std::string(w[1])] = std::move(*p);
            cell_order.emplace_back(w[1]);
        } else if (m.mesh == MeshKind::extruded && e.key == "sweep") {
            m.swept = on_off(e.value).value_or(false);
        }
    }
    if (m.mesh == MeshKind::layered && !domain) return refuse("the layered deck has no domain");
    if (m.mesh == MeshKind::extruded) {
        if (stations.size() < 2 || !std::ranges::is_sorted(stations) || cells.empty()) return refuse("an extrusion the modeller cannot read");
        Pt lo = cells.begin()->second.front(), hi = lo;
        for (const auto& [name, p] : cells) {
            for (const Pt& q : p) {
                lo = {std::min(lo.x, q.x), std::min(lo.y, q.y)};
                hi = {std::max(hi.x, q.x), std::max(hi.y, q.y)};
            }
        }
        m.domain = normalized(from_uv(axis, lo, stations.front()), from_uv(axis, hi, stations.back()));
        // A layered deck fills its domain: the cells must tile their bounds (conforming cells, so by area) and every
        // section of every cell hold a region.
        std::int64_t area = 0;
        for (const auto& [name, p] : cells) area += std::abs(twice_area(p));
        if (area != 2 * (hi.x - lo.x) * (hi.y - lo.y)) m.convert_why = "the extrusion's cells do not fill their bounds";
    }

    std::map<std::string, std::vector<bool>, std::less<>> assigned;  // extruded: [cell][section]
    int region_index = 0;
    for (const auto& s : d->sections) {
        if (s.kind != "region") continue;
        Region r{s.name, {}, {}, region_index++};
        for (const auto& e : s.entries) {
            if (e.key == "material") r.material = std::string(trim(e.value));
            if (e.key == "overlap") r.fill = trim(e.value) == "fill";
            if (m.mesh == MeshKind::layered && is_solid_key(e.key)) {
                auto solid = parse_solid(e, params);
                if (!solid) return refuse(at_line(e, "a solid the modeller cannot read"));
                r.solids.push_back(std::move(*solid));
            } else if (m.mesh == MeshKind::extruded && e.key == "extrude") {
                // extrude = <cells> <t_lo>:<t_hi>: a prism per cell.
                const auto w = words(e.value);
                const std::size_t colon = w.empty() ? std::string_view::npos : w.back().find(':');
                if (w.size() < 2 || colon == std::string_view::npos) return refuse(at_line(e, "an extrude line the modeller cannot read"));
                const auto t0 = parse_length(w.back().substr(0, colon)), t1 = parse_length(w.back().substr(colon + 1));
                if (!t0 || !t1) return refuse(at_line(e, "an extrude line the modeller cannot read"));
                const auto s0 = std::ranges::find(stations, *t0), s1 = std::ranges::find(stations, *t1);
                for (std::size_t k = 0; k + 1 < w.size(); ++k) {
                    const auto c = cells.find(w[k]);
                    if (c == cells.end()) return refuse(at_line(e, "no cell named '" + std::string(w[k]) + "'"));
                    Solid p;
                    p.kind = SolidKind::prism;
                    p.axis = axis;
                    p.from = *t0;
                    p.to = *t1;
                    p.polygon = c->second;
                    r.solids.push_back(std::move(p));
                    auto& row = assigned[c->first];
                    row.resize(stations.size() - 1, false);
                    if (s0 == stations.end() || s1 == stations.end()) {
                        m.convert_why = at_line(e, "a range whose ends the modeller cannot match to stations");
                        continue;
                    }
                    for (auto j = s0 - stations.begin(); j < s1 - stations.begin(); ++j) row[static_cast<std::size_t>(j)] = true;
                }
            }
        }
        if (m.mesh == MeshKind::extruded && r.solids.empty() && m.convert_why.empty()) {
            m.convert_why = "region '" + r.name + "' holds no cell (a layered deck takes one region without solids, its background)";
        }
        m.regions.push_back(std::move(r));
    }
    if (m.mesh == MeshKind::extruded && m.convert_why.empty()) {
        for (const auto& name : cell_order) {
            auto row = assigned[name];
            row.resize(stations.size() - 1, false);
            if (const auto j = std::ranges::find(row, false); j != row.end()) {
                m.convert_why = "cell '" + name + "' is void in section " + std::to_string(j - row.begin()) +
                                " (a layered deck fills its domain)";
                break;
            }
        }
    }

    int contact_index = 0;
    for (const auto& s : d->sections) {
        if (s.kind != "contact") continue;
        Contact c{s.name, {}, {}, true, contact_index++};
        if (const auto* k = find_entry(s, "kind")) c.kind = std::string(trim(k->value));
        const deck::Entry* patch = find_entry(s, "patch");
        const auto it = patch ? patches.find(trim(patch->value)) : patches.end();
        if (it != patches.end()) c.box = it->second;
        else c.editable = false;  // a face's name, or no patch
        m.contacts.push_back(std::move(c));
    }
    return m;
}

// --- writing -------------------------------------------------------------------------------------------------------

std::expected<std::string, std::string> write_model(std::string_view text, const Model& m) {
    if (m.mesh == MeshKind::unsupported) return std::unexpected("this deck's geometry is not editable: " + m.why);
    if (m.mesh == MeshKind::extruded) return std::string(text);
    const auto d = deck::parse(text);
    if (!d) return std::unexpected(d.error().message);
    Rewrite w(text);
    // Unit 26m MT8: a line kept when, read with the model's parameters, it gives the model's solid.
    const Params mp = deck_params(m);
    // [parameters]: a changed expression replaced in place, a new parameter added after the last; a new section
    // before the first.
    const auto param_line = [](const Param& p) { return p.name + " = " + p.expr; };
    const auto ps = std::ranges::find(d->sections, std::string_view("parameters"), &deck::Section::kind);
    if (ps != d->sections.end()) {
        SectionEdit e;
        for (const Param& p : m.parameters) {
            const auto at = std::ranges::find(ps->entries, p.name, &deck::Entry::key);
            if (at == ps->entries.end()) {
                e.after.emplace(ps->entries.empty() ? -1 : static_cast<int>(ps->entries.size()) - 1, param_line(p));
            } else if (std::string(trim(at->value)) != p.expr) {
                // Unit 26n NQ1: the line's label comment kept.
                std::string_view rest = text;
                for (int k = 1; k < at->line; ++k) rest.remove_prefix(std::min(rest.size(), rest.find('\n') + 1));
                std::string_view raw = rest.substr(0, rest.find('\n'));
                if (raw.ends_with('\r')) raw.remove_suffix(1);
                e.replace[static_cast<std::size_t>(at - ps->entries.begin())] = dims::with_value(raw, p.expr);
            }
        }
        w.change(*ps, w.rewrite(*ps, e));
    } else if (!m.parameters.empty()) {
        std::vector<std::string> lines{"[parameters]"};
        for (const Param& p : m.parameters) lines.push_back(param_line(p));
        lines.emplace_back("");
        w.insert(d->sections.empty() ? w.size() : d->sections.front().line - 1, std::move(lines));
    }
    std::vector<const deck::Section*> regions, contacts;
    const deck::Section* mesh = nullptr;
    for (const auto& s : d->sections) {
        if (s.kind == "region") regions.push_back(&s);
        if (s.kind == "contact") contacts.push_back(&s);
        if (s.kind == "mesh") mesh = &s;
    }
    if (mesh == nullptr) return std::unexpected("the deck has no [mesh] section");
    const deck::Entry* type = find_entry(*mesh, "type");
    const bool from_extrusion = type && trim(type->value) == "extrude3d";
    bool swept = false;
    if (const deck::Entry* e = find_entry(*mesh, "sweep")) swept = from_extrusion && on_off(e->value).value_or(false);

    // Regions: slot k (the k-th region section's place) takes the k-th region. A converted extrusion's extrude lines
    // are its solid lines (each replaced by the next prism, the rest after the last).
    const auto solid_entry = [&](std::string_view k) { return is_solid_key(k) || (from_extrusion && k == "extrude"); };
    const auto region_block = [&](const Region& r) -> std::vector<std::string> {
        if (r.origin < 0 || r.origin >= static_cast<int>(regions.size())) {
            std::vector<std::string> out{"", "[region " + r.name + "]"};
            if (!r.material.empty()) out.push_back("material = " + r.material);
            if (r.fill) out.emplace_back("overlap = fill");
            for (const Solid& s : r.solids) out.push_back(solid_line(s, mp));
            return out;
        }
        const deck::Section& s = *regions[static_cast<std::size_t>(r.origin)];
        SectionEdit e;
        if (r.name != s.name) e.header = "[region " + r.name + "]";
        std::size_t solid = 0;
        int last_solid = -1;
        const int last_entry = s.entries.empty() ? -1 : static_cast<int>(s.entries.size()) - 1;
        bool material = false, overlap = false;
        int material_at = -1;
        for (std::size_t k = 0; k < s.entries.size(); ++k) {
            const auto& en = s.entries[k];
            if (en.key == "material") {
                material = true;
                material_at = static_cast<int>(k);
                if (std::string(trim(en.value)) != r.material) {
                    e.replace[k] = r.material.empty() ? std::nullopt : std::optional<std::string>("material = " + r.material);
                }
            } else if (en.key == "overlap") {
                overlap = true;
                const bool says_fill = trim(en.value) == "fill";
                if (says_fill != r.fill) e.replace[k] = r.fill ? std::optional<std::string>("overlap = fill") : std::nullopt;
            } else if (solid_entry(en.key)) {
                last_solid = static_cast<int>(k);
                if (solid >= r.solids.size()) {
                    e.replace[k] = std::nullopt;
                } else if (en.key == "extrude" || parse_solid(en, mp) != r.solids[solid]) {
                    e.replace[k] = solid_line(r.solids[solid], mp);
                }
                ++solid;
            }
        }
        if (!material && !r.material.empty()) e.after.emplace(-1, "material = " + r.material);
        if (!overlap && r.fill) e.after.emplace(material ? material_at : -1, "overlap = fill");
        for (; solid < r.solids.size(); ++solid) e.after.emplace(last_solid >= 0 ? last_solid : last_entry, solid_line(r.solids[solid], mp));
        return w.rewrite(s, e);
    };
    for (std::size_t k = 0; k < std::max(regions.size(), m.regions.size()); ++k) {
        if (k < regions.size() && k < m.regions.size()) {
            auto block = region_block(m.regions[k]);
            if (m.regions[k].origin < 0) block.erase(block.begin());  // a slot keeps its own spacing
            w.change(*regions[k], std::move(block));
        } else if (k < regions.size()) {
            w.drop(*regions[k]);
        } else {
            const int at = regions.empty() ? section_last(w.lines(), *mesh) + 1 : section_last(w.lines(), *regions.back()) + 1;
            Region fresh = m.regions[k];
            fresh.origin = -1;
            w.insert(at, region_block(fresh));
        }
    }

    // Each contact's patch name: a kept contact's as its section says; a new one's its name, made unique among the
    // deck's patches.
    std::vector<std::string> patch_of(m.contacts.size());
    std::set<std::string, std::less<>> taken, deck_used;
    for (const auto& en : mesh->entries) {
        if (en.key.starts_with("patch ")) {
            const auto ws = words(en.key);
            if (ws.size() == 2) taken.emplace(ws[1]);
        }
    }
    for (const deck::Section* s : contacts) {
        if (const auto* p = find_entry(*s, "patch")) deck_used.emplace(trim(p->value));
    }
    for (std::size_t i = 0; i < m.contacts.size(); ++i) {
        const Contact& c = m.contacts[i];
        if (c.origin >= 0 && c.origin < static_cast<int>(contacts.size())) {
            if (const auto* p = find_entry(*contacts[static_cast<std::size_t>(c.origin)], "patch")) patch_of[i] = std::string(trim(p->value));
            taken.insert(patch_of[i]);
        }
    }
    for (std::size_t i = 0; i < m.contacts.size(); ++i) {
        if (!patch_of[i].empty()) continue;
        std::string name = m.contacts[i].name;
        for (int n = 2; taken.contains(name); ++n) name = m.contacts[i].name + "_" + std::to_string(n);
        patch_of[i] = name;
        taken.insert(name);
    }

    // [mesh]: the type (a converted extrusion), layered, domain, the patch boxes.
    {
        SectionEdit e;
        int type_entry = -1, domain_entry = -1, last_patch = -1;
        bool layered = false;
        std::map<std::string, std::size_t, std::less<>> patch_entries;
        static const std::set<std::string_view> extrusion_keys{"axis", "stations", "sweep"};
        static const std::set<std::string_view> swept_keys{"lines", "straddle", "layer", "min_angle", "box"};
        for (std::size_t k = 0; k < mesh->entries.size(); ++k) {
            const auto& en = mesh->entries[k];
            if (en.key == "type") {
                type_entry = static_cast<int>(k);
                if (from_extrusion) e.replace[k] = "type = generate3d";
            } else if (en.key == "layered") {
                layered = true;
                if (!on_off(en.value).value_or(false)) e.replace[k] = "layered = on";
            } else if (en.key == "domain") {
                domain_entry = static_cast<int>(k);
            } else if (en.key.starts_with("patch ")) {
                const auto ws = words(en.key);
                if (ws.size() == 2) patch_entries[std::string(ws[1])] = k;
                last_patch = static_cast<int>(k);
            } else if (from_extrusion && (extrusion_keys.contains(en.key) || en.key.starts_with("cell ") || (swept && swept_keys.contains(en.key)))) {
                e.replace[k] = std::nullopt;
            }
        }
        if (!layered) e.after.emplace(type_entry, "layered = on");
        const std::string domain_line = "domain = " + box_text(m.domain, m.domain_refs, mp);
        if (domain_entry >= 0) {
            std::vector<std::string> was_refs;
            const auto was = parse_box3(mesh->entries[static_cast<std::size_t>(domain_entry)].value, mp, was_refs);
            if (was != m.domain || was_refs != m.domain_refs) e.replace[static_cast<std::size_t>(domain_entry)] = domain_line;
        } else {
            e.after.emplace(type_entry, domain_line);
        }
        // Each contact's patch box; one a removed contact used and no contact uses now is dropped; a new one added.
        std::set<std::string, std::less<>> used;
        for (std::size_t i = 0; i < m.contacts.size(); ++i) {
            const Contact& c = m.contacts[i];
            const std::string& pname = patch_of[i];
            if (!c.editable) {
                used.insert(pname);
                continue;
            }
            const std::string line = "patch " + pname + " = " + patch_text(c.box);
            if (const auto it = patch_entries.find(pname); it != patch_entries.end()) {
                if (!used.contains(pname) && parse_patch(mesh->entries[it->second].value, mp) != c.box) e.replace[it->second] = line;
            } else {
                e.after.emplace(last_patch >= 0 ? last_patch : type_entry, line);
            }
            used.insert(pname);
        }
        for (const auto& [name, k] : patch_entries) {
            if (deck_used.contains(name) && !used.contains(name)) e.replace[k] = std::nullopt;
        }
        w.change(*mesh, w.rewrite(*mesh, e));
    }

    // Contacts: kept ones renamed or given a kind in place; removed ones dropped; new ones appended.
    std::vector<bool> kept(contacts.size(), false);
    std::vector<std::string> appended;
    for (std::size_t i = 0; i < m.contacts.size(); ++i) {
        const Contact& c = m.contacts[i];
        if (c.origin < 0 || c.origin >= static_cast<int>(contacts.size())) {
            appended.insert(appended.end(), {"", "[contact " + c.name + "]", "kind = " + (c.kind.empty() ? std::string("ohmic") : c.kind),
                                             "patch = " + patch_of[i]});
            continue;
        }
        kept[static_cast<std::size_t>(c.origin)] = true;
        const deck::Section& s = *contacts[static_cast<std::size_t>(c.origin)];
        SectionEdit e;
        if (c.name != s.name) e.header = "[contact " + c.name + "]";
        for (std::size_t k = 0; k < s.entries.size(); ++k) {
            if (s.entries[k].key == "kind" && std::string(trim(s.entries[k].value)) != c.kind && !c.kind.empty()) e.replace[k] = "kind = " + c.kind;
        }
        w.change(s, w.rewrite(s, e));
    }
    for (std::size_t k = 0; k < contacts.size(); ++k) {
        if (!kept[k]) w.drop(*contacts[k]);
    }
    if (!appended.empty()) w.insert(w.size(), std::move(appended));
    return w.apply();
}

std::expected<std::string, std::string> convert_to_layered(Model& m) {
    if (m.mesh != MeshKind::extruded) return std::unexpected(std::string("only an extrude3d deck is converted"));
    if (!m.convert_why.empty()) return std::unexpected(m.convert_why);
    m.mesh = MeshKind::layered;
    std::string note = "the mesh becomes a layered generate3d mesh over the domain " + box_text(m.domain) +
                       "; each region's extrude lines become prisms; the extrusion's axis, stations and cells are dropped";
    if (m.swept) note += "; the swept mesh (sweep = on) and its keys are dropped: the polyhedral generator meshes it";
    m.swept = false;
    // A contact on an end plane whose box cuts across a cell's end (slanted cells mostly): the extrusion named only the
    // cells inside the box, a layered deck names the part of every face inside it (slice 4).
    for (const Contact& c : m.contacts) {
        for (const Region& r : m.regions) {
            const bool cut = std::ranges::any_of(r.solids, [&](const Solid& s) {
                if (s.kind != SolidKind::prism || c.box.lo[s.axis] != c.box.hi[s.axis]) return false;
                const Coord plane = c.box.lo[s.axis];
                if (plane != s.from && plane != s.to) return false;
                const Pt lo = uv_of(s.axis, c.box.lo), hi = uv_of(s.axis, c.box.hi);
                Pt plo = s.polygon[0], phi = s.polygon[0];
                for (const Pt& p : s.polygon) {
                    plo = {std::min(plo.x, p.x), std::min(plo.y, p.y)};
                    phi = {std::max(phi.x, p.x), std::max(phi.y, p.y)};
                }
                const bool meets = plo.x < hi.x && phi.x > lo.x && plo.y < hi.y && phi.y > lo.y;
                const bool within = plo.x >= lo.x && phi.x <= hi.x && plo.y >= lo.y && phi.y <= hi.y;
                return meets && !within;
            });
            if (cut) {
                note += "; contact '" + c.name + "' cuts across the end of region '" + r.name +
                        "': a layered deck names the part of a face inside a contact's box, the extrusion only the cells inside it";
            }
        }
    }
    return note;
}

mesh::LayeredGeometry3D to_layered(const Model& m) {
    const auto cm = [](Coord q) { return static_cast<double>(q) * geo::quantum_cm; };
    const auto point = [&](P3 p) { return mesh::Point{cm(p.x), cm(p.y), cm(p.z)}; };
    mesh::LayeredGeometry3D l;
    l.lo = point(m.domain.lo);
    l.hi = point(m.domain.hi);
    l.background = m.background();
    for (std::size_t r = 0; r < m.regions.size(); ++r) {
        for (const Solid& s : m.regions[r].solids) {
            mesh::LayeredSolid v;
            v.region = static_cast<std::int32_t>(r);
            v.fill = m.regions[r].fill;
            if (s.kind == SolidKind::box) {
                v.kind = mesh::LayeredSolid::Kind::box;
                v.lo = point(s.box.lo);
                v.hi = point(s.box.hi);
            } else {
                v.kind = mesh::LayeredSolid::Kind::prism;
                v.axis = s.axis;
                v.from = cm(s.from);
                v.to = cm(s.to);
                for (const Pt& p : s.polygon) v.polygon.push_back({cm(p.x), cm(p.y), 0.0});
            }
            l.solids.push_back(std::move(v));
        }
    }
    for (const Contact& c : m.contacts) {
        if (c.editable) l.patches.push_back({c.name, point(c.box.lo), point(c.box.hi)});
    }
    return l;
}

// --- edits ---------------------------------------------------------------------------------------------------------

namespace {

bool valid(const Model& m, SolidRef s) {
    return m.editable() && s.region >= 0 && s.region < static_cast<int>(m.regions.size()) && s.solid >= 0 &&
           s.solid < static_cast<int>(m.regions[static_cast<std::size_t>(s.region)].solids.size());
}

Solid& at(Model& m, SolidRef s) { return m.regions[static_cast<std::size_t>(s.region)].solids[static_cast<std::size_t>(s.solid)]; }

bool valid_region(const Model& m, int region) { return m.editable() && region >= 0 && region < static_cast<int>(m.regions.size()); }

std::string fresh(const Model& m, std::string base) {
    if (!identifier(base)) base = "region";
    for (int n = 1;; ++n) {
        std::string name = base + "_" + std::to_string(n);
        if (std::ranges::none_of(m.regions, [&](const Region& r) { return r.name == name; })) return name;
    }
}

int add_region(Model& m, Solid s, std::string material) {
    Region r{fresh(m, material), std::move(material), {std::move(s)}, -1};
    m.regions.push_back(std::move(r));
    return static_cast<int>(m.regions.size()) - 1;
}

bool solid_ok(const Solid& s) {
    if (s.kind == SolidKind::box) return box_ok(s.box);
    return s.axis >= 0 && s.axis <= 2 && s.from < s.to && within(s.from) && within(s.to) && prism_ok(s.polygon);
}

// `s` replaced by `t` when `t` is valid.
bool replace(Solid& s, Solid t) {
    if (!solid_ok(t)) return false;
    s = std::move(t);
    return true;
}

}  // namespace

bool prism_ok(const std::vector<Pt>& p) noexcept { return geo::polygon_ok(p); }  // slanted edges too (slice 4, SD6)

int add_box(Model& m, Box3 b, std::string material) {
    if (!m.editable() || !box_ok(b)) return -1;
    return add_region(m, Solid{SolidKind::box, b, 2, 0, 0, {}}, std::move(material));
}

int add_prism(Model& m, int axis, Coord from, Coord to, std::vector<Pt> polygon, std::string material) {
    Solid s{SolidKind::prism, {}, axis, from, to, std::move(polygon)};
    if (!m.editable() || !solid_ok(s)) return -1;
    return add_region(m, std::move(s), std::move(material));
}

bool move_solid(Model& m, SolidRef s, P3 delta) {
    if (!valid(m, s)) return false;
    Solid t = at(m, s);
    if (t.kind == SolidKind::box) {
        for (int a = 0; a < 3; ++a) {
            t.box.lo[a] += delta[a];
            t.box.hi[a] += delta[a];
        }
    } else {
        t.from += delta[t.axis];
        t.to += delta[t.axis];
        const Pt d = uv_of(t.axis, delta);
        for (Pt& p : t.polygon) p = {p.x + d.x, p.y + d.y};
    }
    return replace(at(m, s), std::move(t));
}

bool set_box(Model& m, SolidRef s, Box3 b) {
    if (!valid(m, s) || at(m, s).kind != SolidKind::box) return false;
    Solid t = at(m, s);
    t.box = b;
    return replace(at(m, s), std::move(t));
}

bool set_extent(Model& m, SolidRef s, Coord from, Coord to) {
    if (!valid(m, s) || at(m, s).kind != SolidKind::prism) return false;
    Solid t = at(m, s);
    t.from = from;
    t.to = to;
    return replace(at(m, s), std::move(t));
}

bool move_face(Model& m, SolidRef s, int face, Coord offset) {
    if (!valid(m, s) || face < 0) return false;
    Solid t = at(m, s);
    if (t.kind == SolidKind::box) {
        if (face >= 6) return false;
        Coord& c = (face % 2 == 0) ? t.box.lo[face / 2] : t.box.hi[face / 2];
        c += offset;
    } else if (face == 0) {
        t.from += offset;
    } else if (face == 1) {
        t.to += offset;
    } else {
        const std::size_t n = t.polygon.size(), k = static_cast<std::size_t>(face - 2);
        if (k >= n) return false;
        Pt& a = t.polygon[k];
        Pt& b = t.polygon[(k + 1) % n];
        if (a.y == b.y) {  // parallel to u: moves along v
            a.y += offset;
            b.y += offset;
        } else if (a.x == b.x) {
            a.x += offset;
            b.x += offset;
        } else {
            return false;  // a slanted side: slice 4
        }
    }
    return replace(at(m, s), std::move(t));
}

bool set_vertex(Model& m, SolidRef s, int k, Pt p, bool follow) {
    if (!valid(m, s) || at(m, s).kind != SolidKind::prism) return false;
    Solid t = at(m, s);
    auto& w = t.polygon;
    if (k < 0 || k >= static_cast<int>(w.size())) return false;
    const std::size_t n = w.size(), i = static_cast<std::size_t>(k);
    const Pt old = w[i];
    if (follow) {
        for (const std::size_t j : {(i + n - 1) % n, (i + 1) % n}) {
            if (w[j].x == old.x) w[j].x = p.x;
            else if (w[j].y == old.y) w[j].y = p.y;
        }
    }
    w[i] = p;
    return replace(at(m, s), std::move(t));
}

bool insert_vertex(Model& m, SolidRef s, int after, Pt p) {
    if (!valid(m, s) || at(m, s).kind != SolidKind::prism) return false;
    Solid t = at(m, s);
    auto& v = t.polygon;
    if (after < 0 || after >= static_cast<int>(v.size())) return false;
    const Pt a = v[static_cast<std::size_t>(after)], b = v[(static_cast<std::size_t>(after) + 1) % v.size()];
    // Exactly on the open edge, slanted edges too (64-bit cross product of quanta within the coordinate limit).
    const std::int64_t cross = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
    const bool on = cross == 0 && p != a && p != b && p.x >= std::min(a.x, b.x) && p.x <= std::max(a.x, b.x) &&
                    p.y >= std::min(a.y, b.y) && p.y <= std::max(a.y, b.y);
    if (!on) return false;
    v.insert(v.begin() + after + 1, p);
    return replace(at(m, s), std::move(t));
}

bool remove_vertex(Model& m, SolidRef s, int k) {
    if (!valid(m, s) || at(m, s).kind != SolidKind::prism) return false;
    Solid t = at(m, s);
    if (k < 0 || k >= static_cast<int>(t.polygon.size())) return false;
    t.polygon.erase(t.polygon.begin() + k);
    return replace(at(m, s), std::move(t));
}

bool delete_solid(Model& m, SolidRef s) {
    if (!valid(m, s)) return false;
    auto& solids = m.regions[static_cast<std::size_t>(s.region)].solids;
    solids.erase(solids.begin() + s.solid);
    if (solids.empty()) m.regions.erase(m.regions.begin() + s.region);
    return true;
}

bool rename_region(Model& m, int region, std::string name) {
    if (!valid_region(m, region) || !identifier(name)) return false;
    for (std::size_t k = 0; k < m.regions.size(); ++k) {
        if (static_cast<int>(k) != region && m.regions[k].name == name) return false;
    }
    m.regions[static_cast<std::size_t>(region)].name = std::move(name);
    return true;
}

bool set_material(Model& m, int region, std::string material) {
    if (!valid_region(m, region) || !identifier(material)) return false;
    m.regions[static_cast<std::size_t>(region)].material = std::move(material);
    return true;
}

bool raise_region(Model& m, int region) {
    if (!valid_region(m, region) || region + 1 >= static_cast<int>(m.regions.size())) return false;
    std::swap(m.regions[static_cast<std::size_t>(region)], m.regions[static_cast<std::size_t>(region + 1)]);
    return true;
}

bool lower_region(Model& m, int region) {
    if (!valid_region(m, region) || region == 0) return false;
    return raise_region(m, region - 1);
}

int duplicate_solid(Model& m, SolidRef s, P3 delta) {
    if (!valid(m, s)) return -1;
    const std::string material = m.regions[static_cast<std::size_t>(s.region)].material;
    const int r = add_region(m, at(m, s), material);
    if (!move_solid(m, {r, 0}, delta)) {
        m.regions.pop_back();
        return -1;
    }
    return r;
}

bool mirror_solid(Model& m, SolidRef s, int axis) {
    if (!valid(m, s) || axis < 0 || axis > 2) return false;
    Solid& t = at(m, s);
    if (t.kind == SolidKind::box || axis == t.axis) return true;  // symmetric about its centre across `axis`
    const Box3 b = t.bounds();
    const bool u = axis == (t.axis + 1) % 3;
    for (Pt& p : t.polygon) {
        if (u) p.x = b.lo[axis] + b.hi[axis] - p.x;
        else p.y = b.lo[axis] + b.hi[axis] - p.y;
    }
    std::ranges::reverse(t.polygon);  // keep the orientation
    return true;
}


// --- Unit 26m slice 1 (MT4, MT5) ---

bool set_overlap(Model& m, int region, bool fill) {
    if (!valid_region(m, region)) return false;
    Region& r = m.regions[static_cast<std::size_t>(region)];
    if (r.solids.empty() || r.fill == fill) return false;
    r.fill = fill;
    return true;
}

namespace {

Coord floor_half(Coord v) noexcept { return v >= 0 ? v / 2 : -((-v + 1) / 2); }

double turn_of(const std::vector<Pt>& p) {
    double a = 0.0;
    for (std::size_t k = 0; k < p.size(); ++k) {
        const Pt& u = p[k];
        const Pt& w = p[(k + 1) % p.size()];
        a += static_cast<double>(u.x) * static_cast<double>(w.y) - static_cast<double>(w.x) * static_cast<double>(u.y);
    }
    return a;
}

// `s` under a map of the quanta's lattice onto itself (axes permuted, signs flipped, shifted); a prism's polygon kept
// in its orientation.
template <class Map>
Solid mapped(const Solid& s, Map&& map) {
    Solid t = s;
    if (s.kind == SolidKind::box) {
        t.box = normalized(map(s.box.lo), map(s.box.hi));
        return t;
    }
    const P3 a0 = map(from_uv(s.axis, s.polygon[0], s.from)), a1 = map(from_uv(s.axis, s.polygon[0], s.to));
    int axis = 0;
    while (axis < 2 && a0[axis] == a1[axis]) ++axis;
    t.axis = axis;
    t.from = std::min(a0[axis], a1[axis]);
    t.to = std::max(a0[axis], a1[axis]);
    t.polygon.clear();
    for (const Pt& uv : s.polygon) t.polygon.push_back(uv_of(axis, map(from_uv(s.axis, uv, s.from))));
    if ((turn_of(t.polygon) > 0) != (turn_of(s.polygon) > 0)) std::ranges::reverse(t.polygon);
    return t;
}

}  // namespace

bool mirror_solid_at(Model& m, SolidRef s, int axis, Coord twice) {
    if (!valid(m, s) || axis < 0 || axis > 2) return false;
    return replace(at(m, s), mapped(at(m, s), [&](P3 p) {
                       p[axis] = twice - p[axis];
                       return p;
                   }));
}

bool rotate_solid(Model& m, SolidRef s, int axis) {
    if (!valid(m, s) || axis < 0 || axis > 2) return false;
    const Box3 b = at(m, s).bounds();
    const int u = (axis + 1) % 3, v = (axis + 2) % 3;
    const Coord cu2 = b.lo[u] + b.hi[u], cv2 = b.lo[v] + b.hi[v];
    const Coord su = floor_half(cu2 + cv2), sv = floor_half(cv2 - cu2);  // u' = cu + cv - v, v' = cv - cu + u
    return replace(at(m, s), mapped(at(m, s), [&](P3 p) {
                       const Coord pu = p[u], pv = p[v];
                       p[u] = su - pv;
                       p[v] = sv + pu;
                       return p;
                   }));
}

int array_solid(Model& m, SolidRef s, P3 pitch, int count) {
    if (!valid(m, s) || count < 2 || count > 100 || pitch == P3{}) return -1;
    const Model before = m;
    const bool fill = m.regions[static_cast<std::size_t>(s.region)].fill;
    int first = -1;
    for (int k = 1; k < count; ++k) {
        const int r = duplicate_solid(m, s, {k * pitch.x, k * pitch.y, k * pitch.z});
        if (r < 0) {
            m = before;
            return -1;
        }
        m.regions[static_cast<std::size_t>(r)].fill = fill;
        if (first < 0) first = r;
    }
    return first;
}

bool align_solids(Model& m, const std::vector<SolidRef>& picked, int face, Coord to) {
    if (face < 0 || face > 5 || picked.empty()) return false;
    std::vector<SolidRef> solids;  // each once
    for (const SolidRef s : picked) {
        if (!valid(m, s)) return false;
        if (std::ranges::find(solids, s) == solids.end()) solids.push_back(s);
    }
    const int axis = face / 2;
    const Model before = m;
    bool moved = false;
    for (const SolidRef s : solids) {
        const Box3 b = at(m, s).bounds();
        P3 d{};
        d[axis] = to - (face % 2 == 0 ? b.lo[axis] : b.hi[axis]);
        if (d == P3{}) continue;
        if (!move_solid(m, s, d)) {
            m = before;
            return false;
        }
        moved = true;
    }
    return moved;
}

// --- Unit 26m slice 2 (MT6, MT7, MT12) ---

std::optional<EdgeCorner> edge_corner(const Solid& s, int edge) {
    const auto list = s.edges();
    if (edge < 0 || edge >= static_cast<int>(list.size())) return std::nullopt;
    if (s.kind == SolidKind::prism) {
        if (edge % 3 != 2) return std::nullopt;  // a cap's edge
        return EdgeCorner{s, edge / 3};
    }
    const auto [i, j] = list[static_cast<std::size_t>(edge)];
    const int bit = i ^ j, axis = bit == 1 ? 0 : bit == 2 ? 1 : 2, u = (axis + 1) % 3, v = (axis + 2) % 3;
    const Box3& b = s.box;
    Solid p{SolidKind::prism, {}, axis, b.lo[axis], b.hi[axis],
            {{b.lo[u], b.lo[v]}, {b.hi[u], b.lo[v]}, {b.hi[u], b.hi[v]}, {b.lo[u], b.hi[v]}}};
    const bool hu = (i >> u) & 1, hv = (i >> v) & 1;
    return EdgeCorner{std::move(p), hv ? (hu ? 2 : 3) : (hu ? 1 : 0)};
}

namespace {

template <class Cut>
bool cut_edge(Model& m, SolidRef s, int edge, Cut&& cut) {
    if (!valid(m, s)) return false;
    auto ec = edge_corner(at(m, s), edge);
    if (!ec) return false;
    auto poly = cut(ec->prism.polygon, ec->corner);
    if (!poly) return false;
    ec->prism.polygon = std::move(*poly);
    return replace(at(m, s), std::move(ec->prism));
}

}  // namespace

bool chamfer_edge(Model& m, SolidRef s, int edge, Coord d) {
    return cut_edge(m, s, edge, [&](const std::vector<Pt>& p, int k) { return geo::chamfer_corner(p, k, d); });
}

bool fillet_edge(Model& m, SolidRef s, int edge, Coord r, int segments) {
    return cut_edge(m, s, edge, [&](const std::vector<Pt>& p, int k) { return geo::fillet_corner(p, k, r, segments); });
}

bool fillet_vertex(Model& m, SolidRef s, int k, Coord r, int segments) {
    if (!valid(m, s) || at(m, s).kind != SolidKind::prism) return false;
    const int n = static_cast<int>(at(m, s).polygon.size());
    if (k < 0 || k >= n) return false;
    return fillet_edge(m, s, 3 * k + 2, r, segments);
}

bool taper_edge(Model& m, SolidRef s, int edge, Pt delta) {
    if (!valid(m, s) || delta == Pt{}) return false;
    auto ec = edge_corner(at(m, s), edge);
    if (!ec) return false;
    Pt& p = ec->prism.polygon[static_cast<std::size_t>(ec->corner)];
    p = {p.x + delta.x, p.y + delta.y};
    return replace(at(m, s), std::move(ec->prism));
}

int add_cylinder(Model& m, int axis, Coord from, Coord to, Pt centre, Coord r, int sides, std::string material) {
    if (r <= 0 || sides < 3 || sides > 360 || axis < 0 || axis > 2) return -1;
    return add_prism(m, axis, std::min(from, to), std::max(from, to), geo::regular_polygon(centre, r, sides), std::move(material));
}

std::optional<Normal> face_normal(const Solid& s, int face) {
    if (s.kind == SolidKind::box) {
        if (face < 0 || face >= 6) return std::nullopt;
        return Normal{face / 2, face % 2 ? 1 : -1};
    }
    if (face == 0 || face == 1) return Normal{s.axis, face == 1 ? 1 : -1};
    const std::size_t n = s.polygon.size(), k = static_cast<std::size_t>(face - 2);
    if (face < 2 || k >= n) return std::nullopt;
    const Pt a = s.polygon[k], b = s.polygon[(k + 1) % n];
    const int o = twice_area(s.polygon) > 0 ? 1 : -1;  // counter-clockwise: the outward normal is (dv, -du)
    if (a.y == b.y) return Normal{(s.axis + 2) % 3, (b.x > a.x ? -1 : 1) * o};
    if (a.x == b.x) return Normal{(s.axis + 1) % 3, (b.y > a.y ? 1 : -1) * o};
    return std::nullopt;
}

std::optional<Solid> face_solid(const Solid& s, int face, Coord margin, Coord depth) {
    const auto n = face_normal(s, face);
    if (!n || margin < 0 || depth == 0) return std::nullopt;
    Box3 f;  // the face, flat across n->axis
    if (s.kind == SolidKind::box) {
        f = s.box;
    } else if (face <= 1) {
        const Box3 b = s.bounds();
        const Pt lo = uv_of(s.axis, b.lo), hi = uv_of(s.axis, b.hi);
        const bool rectangle = s.polygon.size() == 4 && std::ranges::all_of(s.polygon, [&](Pt p) {
            return (p.x == lo.x || p.x == hi.x) && (p.y == lo.y || p.y == hi.y);
        });
        if (!rectangle) {
            if (margin != 0) return std::nullopt;
            const Coord c = face == 0 ? s.from : s.to, e = c + n->sign * depth;
            Solid t{SolidKind::prism, {}, s.axis, std::min(c, e), std::max(c, e), s.polygon};
            if (!solid_ok(t)) return std::nullopt;
            return t;
        }
        f = b;
    } else {
        const std::size_t k = static_cast<std::size_t>(face - 2);
        const Pt a = s.polygon[k], b = s.polygon[(k + 1) % s.polygon.size()];
        f = normalized(from_uv(s.axis, a, s.from), from_uv(s.axis, b, s.to));
    }
    const int ax = n->axis;
    const Coord c = n->sign > 0 ? f.hi[ax] : f.lo[ax], e = c + n->sign * depth;
    f.lo[ax] = std::min(c, e);
    f.hi[ax] = std::max(c, e);
    for (const int o : {(ax + 1) % 3, (ax + 2) % 3}) {
        f.lo[o] += margin;
        f.hi[o] -= margin;
    }
    if (!box_ok(f)) return std::nullopt;
    return Solid{SolidKind::box, f, 2, 0, 0, {}};
}

int inset_face(Model& m, SolidRef s, int face, Coord margin, Coord depth, std::string material) {
    if (!valid(m, s)) return -1;
    auto t = face_solid(at(m, s), face, margin, depth);
    if (!t) return -1;
    Region& own = m.regions[static_cast<std::size_t>(s.region)];
    if (material == own.material) {  // its own material: joins its region, so no interface forms between them
        own.solids.push_back(std::move(*t));
        return s.region;
    }
    return add_region(m, std::move(*t), std::move(material));
}

int extrude_face(Model& m, SolidRef s, int face, Coord depth, std::string material) {
    if (depth <= 0) return -1;
    return inset_face(m, s, face, 0, depth, std::move(material));
}

std::expected<double, std::string> automatic_resolution(const Model& m, double max_size) {
    const auto complex = mesh::layered_partition3d(to_layered(m));
    if (!complex) return std::unexpected(complex.error().message);
    mesh::MeshOptions3D o;
    o.max_size = max_size;
    const auto s = mesh::largest_meshable_resolution(*complex, o);
    if (!s) return std::unexpected(s.error().message);
    return *s;
}

std::optional<Extent> extent(const Model& m) {
    const auto complex = mesh::layered_partition3d(to_layered(m));
    if (!complex) return std::nullopt;
    const auto& g = *complex;
    std::vector<std::vector<int>> owners(g.facets.size());  // the regions of the volumes each facet bounds
    std::map<std::pair<std::size_t, std::size_t>, std::vector<std::pair<std::array<double, 3>, std::pair<int, int>>>> edges;
    Extent e;
    for (const auto& v : g.volumes) {
        for (const std::size_t f : v.facets) owners[f].push_back(v.region);
    }
    for (std::size_t f = 0; f < g.facets.size(); ++f) {
        const auto& vs = g.facets[f].vertices;
        std::array<double, 3> n{};  // Newell's normal: twice the area along it
        for (std::size_t i = 0; i < vs.size(); ++i) {
            const auto& a = g.vertices[vs[i]];
            const auto& b = g.vertices[vs[(i + 1) % vs.size()]];
            n[0] += (a[1] - b[1]) * (a[2] + b[2]);
            n[1] += (a[2] - b[2]) * (a[0] + b[0]);
            n[2] += (a[0] - b[0]) * (a[1] + b[1]);
        }
        const double area = 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        const bool interface = owners[f].size() == 2 && owners[f][0] != owners[f][1];
        if (interface) e.interface_area += area;
        if (!(area > 0.0) || (!interface && owners[f].size() != 1)) continue;
        for (double& c : n) c /= 2.0 * area;
        for (std::size_t i = 0; i < vs.size(); ++i) {
            const std::size_t a = vs[i], b = vs[(i + 1) % vs.size()];
            auto& on = edges[{std::min(a, b), std::max(a, b)}];
            on.push_back({n, interface ? std::pair<int, int>(std::minmax(owners[f][0], owners[f][1])) : std::pair<int, int>{-1, -1}});
        }
    }
    // A feature line: an edge of an interface where its facets turn (a bend, or the boundary) or two region pairs meet.
    for (const auto& [ab, on] : edges) {
        const bool any = std::ranges::any_of(on, [](const auto& o) { return o.second.first >= 0; });
        bool feature = false;
        for (std::size_t i = 0; any && i < on.size() && !feature; ++i) {
            for (std::size_t j = i + 1; j < on.size() && !feature; ++j) {
                const auto& p = on[i].first;
                const auto& q = on[j].first;
                const double cx = p[1] * q[2] - p[2] * q[1], cy = p[2] * q[0] - p[0] * q[2], cz = p[0] * q[1] - p[1] * q[0];
                feature = cx * cx + cy * cy + cz * cz > 1e-12 || (on[i].second != on[j].second && on[i].second.first >= 0 && on[j].second.first >= 0);
            }
        }
        if (!feature) continue;
        const auto& a = g.vertices[ab.first];
        const auto& b = g.vertices[ab.second];
        e.feature_length += std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
    }
    const auto report = mesh::check_complex(g);
    if (!report) return std::nullopt;
    for (const double v : report->volume) e.volume += v;
    return e;
}

double node_estimate(const Extent& e, double h, double s) {
    // Least squares on the slice 2 table (7 meshes of tool_deck layouts, 7k-15k nodes): within ±10 % there (in-sample).
    constexpr double a = 2.65, b = 1.12, c = 13.3;
    if (!(h > 0.0)) return 0.0;
    const double hp = s > 0.0 ? std::min(h, 4.0 * s) : h;
    return a * e.volume / (h * h * h) + b * e.interface_area / (hp * hp) + c * e.feature_length / hp;
}

std::optional<double> mesh_length(std::string_view text, std::string_view key) {
    const auto d = deck::parse(text);
    if (!d) return std::nullopt;
    const deck::Section* mesh = nullptr;  // the last, as read_model, write_model and set_mesh_length take it
    for (const auto& s : d->sections) {
        if (s.kind == "mesh") mesh = &s;
    }
    if (!mesh) return std::nullopt;
    const deck::Entry* e = find_entry(*mesh, key);
    if (!e) return std::nullopt;
    const auto q = geo::parse_length(trim(e->value));
    if (!q) return std::nullopt;
    return static_cast<double>(*q) / (per_um * 1e4);
}

std::expected<std::string, std::string> set_mesh_length(std::string_view text, std::string_view key, Coord q) {
    const auto d = deck::parse(text);
    if (!d) return std::unexpected(d.error().message);
    const deck::Section* mesh = nullptr;
    for (const auto& s : d->sections) {
        if (s.kind == "mesh") mesh = &s;
    }
    if (!mesh) return std::unexpected(std::string("the deck has no [mesh] section"));
    Lines l = split(text);
    const std::string line = std::string(key) + " = " + geo::format_length(q);
    const deck::Entry* e = find_entry(*mesh, key);
    int at = 0;
    if (e) {
        const int first = e->line - 1, last = entry_last(l, first);
        l.eol[static_cast<std::size_t>(first)] = l.eol[static_cast<std::size_t>(last)];  // the entry's own end (none at EOF)
        l.body.erase(l.body.begin() + first + 1, l.body.begin() + last + 1);
        l.eol.erase(l.eol.begin() + first + 1, l.eol.begin() + last + 1);
        l.body[static_cast<std::size_t>(first)] = line;
    } else {
        at = section_last(l, *mesh) + 1;
        l.body.insert(l.body.begin() + at, line);
        l.eol.insert(l.eol.begin() + at, l.nl);
        if (at > 0 && l.eol[static_cast<std::size_t>(at - 1)].empty()) l.eol[static_cast<std::size_t>(at - 1)] = l.nl;
    }
    std::string out;
    for (std::size_t k = 0; k < l.body.size(); ++k) out += l.body[k] + l.eol[k];
    return out;
}

Coord keep_clear(const Model& m, const std::vector<Solid>& moving, const std::vector<SolidRef>& skip, int axis, Coord delta,
                 Coord margin) {
    if (margin <= 0 || axis < 0 || axis > 2 || moving.empty()) return delta;
    std::vector<std::pair<Coord, Coord>> pairs;  // (a moving face's coordinate, a plane it must keep clear of)
    for (const Solid& mv : moving) {
        const Box3 mb = mv.bounds();
        std::vector<Coord> planes{m.domain.lo[axis], m.domain.hi[axis]};
        for (std::size_t r = 0; r < m.regions.size(); ++r) {
            for (std::size_t k = 0; k < m.regions[r].solids.size(); ++k) {
                if (std::ranges::find(skip, SolidRef{static_cast<int>(r), static_cast<int>(k)}) != skip.end()) continue;
                const Box3 b = m.regions[r].solids[k].bounds();
                bool near = true;
                for (const int o : {(axis + 1) % 3, (axis + 2) % 3}) near = near && b.lo[o] - margin <= mb.hi[o] && mb.lo[o] <= b.hi[o] + margin;
                if (near) {
                    planes.push_back(b.lo[axis]);
                    planes.push_back(b.hi[axis]);
                }
            }
        }
        for (const Coord f : {mb.lo[axis], mb.hi[axis]}) {
            for (const Coord p : planes) pairs.emplace_back(f, p);
        }
    }
    const auto clear = [&](Coord d) {
        return std::ranges::none_of(pairs, [&](const auto& fp) {
            const Coord gap = std::abs(fp.first + d - fp.second);
            return gap > 0 && gap < margin;
        });
    };
    if (clear(delta)) return delta;
    std::optional<Coord> best;
    for (const auto& [f, p] : pairs) {
        for (const Coord d : {p - f, p - f - margin, p - f + margin}) {
            if (std::abs(d - delta) > 2 * margin || !clear(d)) continue;
            if (!best || std::abs(d - delta) < std::abs(*best - delta) || (std::abs(d - delta) == std::abs(*best - delta) && d < *best)) best = d;
        }
    }
    return best.value_or(delta);
}

// --- Unit 26n slice 4 (DT8) ------------------------------------------------------------------------------------------

std::string outside_domain(const Model& m) {
    static constexpr const char* face[3][2] = {{"x_min", "x_max"}, {"y_min", "y_max"}, {"z_min", "z_max"}};
    for (const Region& r : m.regions) {
        for (const Solid& s : r.solids) {
            const Box3 b = s.bounds();
            for (int a = 0; a < 3; ++a) {
                const Coord below = m.domain.lo[a] - b.lo[a], above = b.hi[a] - m.domain.hi[a];
                if (below > 0 || above > 0) {
                    const bool high = above >= below;
                    return "region '" + r.name + "' reaches " + format_length(high ? above : below) + " past the domain's " + face[a][high ? 1 : 0] +
                           " face";
                }
            }
        }
    }
    return {};
}

P3 clamp_to(const Box3& d, P3 p) noexcept {
    for (int a = 0; a < 3; ++a) p[a] = std::clamp(p[a], d.lo[a], d.hi[a]);
    return p;
}

P3 clamp_offset(const Box3& d, const std::vector<Solid>& moving, P3 delta) noexcept {
    if (moving.empty()) return delta;
    for (int a = 0; a < 3; ++a) {
        Coord lo = std::numeric_limits<Coord>::max(), hi = std::numeric_limits<Coord>::lowest();
        for (const Solid& s : moving) {
            const Box3 b = s.bounds();
            lo = std::min(lo, b.lo[a]);
            hi = std::max(hi, b.hi[a]);
        }
        const Coord down = std::min<Coord>(0, d.lo[a] - lo), up = std::max<Coord>(0, d.hi[a] - hi);
        delta[a] = std::clamp(delta[a], down, up);
    }
    return delta;
}

Coord clamp_scalar(Coord v, const std::function<bool(Coord)>& inside) {
    if (v == 0 || inside(v)) return v;
    Coord good = 0, bad = v;  // inside(good), !inside(bad)
    while (bad - good > 1 || good - bad > 1) {
        const Coord mid = good + (bad - good) / 2;
        (inside(mid) ? good : bad) = mid;
    }
    return good;
}

bool selects_boundary(const Box3& d, const Box3& b) noexcept {
    for (int a = 0; a < 3; ++a) {
        if (b.lo[a] > b.hi[a] || b.lo[a] < d.lo[a] || b.hi[a] > d.hi[a]) return false;
    }
    for (int a = 0; a < 3; ++a) {
        const int u = (a + 1) % 3, v = (a + 2) % 3;
        if ((b.lo[a] == d.lo[a] || b.hi[a] == d.hi[a]) && b.lo[u] < b.hi[u] && b.lo[v] < b.hi[v]) return true;
    }
    return false;
}

int add_contact(Model& m, std::string name, std::string kind, Box3 b) {
    if (!m.editable() || !identifier(name) || !selects_boundary(m.domain, b)) return -1;
    if (std::ranges::any_of(m.contacts, [&](const Contact& c) { return c.name == name; })) return -1;
    m.contacts.push_back({std::move(name), std::move(kind), b, true, -1});
    return static_cast<int>(m.contacts.size()) - 1;
}

bool set_contact(Model& m, int contact, Box3 b) {
    if (!m.editable() || contact < 0 || contact >= static_cast<int>(m.contacts.size())) return false;
    Contact& c = m.contacts[static_cast<std::size_t>(contact)];
    if (!c.editable || !selects_boundary(m.domain, b)) return false;
    c.box = b;
    return true;
}

bool move_contact(Model& m, int contact, P3 delta) {
    if (!m.editable() || contact < 0 || contact >= static_cast<int>(m.contacts.size())) return false;
    Box3 b = m.contacts[static_cast<std::size_t>(contact)].box;
    for (int a = 0; a < 3; ++a) {
        b.lo[a] += delta[a];
        b.hi[a] += delta[a];
    }
    return set_contact(m, contact, b);
}

bool delete_contact(Model& m, int contact) {
    if (!m.editable() || contact < 0 || contact >= static_cast<int>(m.contacts.size())) return false;
    m.contacts.erase(m.contacts.begin() + contact);
    return true;
}

bool rename_contact(Model& m, int contact, std::string name) {
    if (!m.editable() || contact < 0 || contact >= static_cast<int>(m.contacts.size()) || !identifier(name)) return false;
    for (std::size_t k = 0; k < m.contacts.size(); ++k) {
        if (static_cast<int>(k) != contact && m.contacts[k].name == name) return false;
    }
    m.contacts[static_cast<std::size_t>(contact)].name = std::move(name);
    return true;
}

bool set_contact_kind(Model& m, int contact, std::string kind) {
    if (!m.editable() || contact < 0 || contact >= static_cast<int>(m.contacts.size()) || !identifier(kind)) return false;
    m.contacts[static_cast<std::size_t>(contact)].kind = std::move(kind);
    return true;
}

bool set_domain(Model& m, Box3 d) {
    if (!m.editable() || !box_ok(d)) return false;
    // A contact's coordinate on a face of the old domain moves with the face; the others are clipped.
    std::vector<Contact> contacts = m.contacts;
    for (Contact& c : contacts) {
        if (!c.editable) continue;
        for (int a = 0; a < 3; ++a) {
            for (Coord* x : {&c.box.lo[a], &c.box.hi[a]}) {
                if (*x == m.domain.lo[a]) *x = d.lo[a];
                else if (*x == m.domain.hi[a]) *x = d.hi[a];
                else *x = std::clamp(*x, d.lo[a], d.hi[a]);
            }
        }
        if (!selects_boundary(d, c.box)) return false;
    }
    m.domain = d;
    m.contacts = std::move(contacts);
    return true;
}

// --- snapping ------------------------------------------------------------------------------------------------------

Coord grid_step(double px_per_um, int dpi) noexcept {
    geo::View v;
    v.px_per_um = px_per_um;
    v.dpi = dpi;
    return geo::grid_step(v);
}

namespace {

double distance(P3 a, P3 b) {
    const double dx = static_cast<double>(a.x - b.x), dy = static_cast<double>(a.y - b.y), dz = static_cast<double>(a.z - b.z);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

Coord round_to(double v) { return static_cast<Coord>(std::llround(v)); }

// The closed polygon contains `p` (boundary included; exact).
bool inside(const std::vector<Pt>& poly, Pt p) {
    bool in = false;
    const std::size_t n = poly.size();
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        const Pt a = poly[i], b = poly[j];
        const std::int64_t cross = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
        if (cross == 0 && p.x >= std::min(a.x, b.x) && p.x <= std::max(a.x, b.x) && p.y >= std::min(a.y, b.y) && p.y <= std::max(a.y, b.y)) {
            return true;
        }
        if ((a.y > p.y) != (b.y > p.y)) {
            // p.x < the crossing's x, compared exactly: (p.x - a.x)(b.y - a.y) < (b.x - a.x)(p.y - a.y), sign by b.y - a.y.
            const std::int64_t l = (p.x - a.x) * (b.y - a.y), r = (b.x - a.x) * (p.y - a.y);
            if ((b.y > a.y) ? l < r : l > r) in = !in;
        }
    }
    return in;
}

}  // namespace

Snap snap(const Model& m, const WorkPlane& plane, P3 raw, Coord tolerance, std::optional<SolidRef> skip, Coord grid) {
    const auto on_grid = [&](Coord c) { return round_to(static_cast<double>(c) / static_cast<double>(grid)) * grid; };
    std::vector<Solid> solids;
    const auto as_box = [](const Box3& b) { return Solid{SolidKind::box, b, 2, 0, 0, {}}; };
    solids.push_back(as_box(m.domain));
    for (std::size_t r = 0; r < m.regions.size(); ++r) {
        for (std::size_t k = 0; k < m.regions[r].solids.size(); ++k) {
            if (skip && *skip == SolidRef{static_cast<int>(r), static_cast<int>(k)}) continue;
            solids.push_back(m.regions[r].solids[k]);
        }
    }
    for (const Contact& c : m.contacts) {
        if (c.editable) solids.push_back(as_box(c.box));
    }
    const double tol = static_cast<double>(tolerance);
    const auto best_of = [&](auto&& each) {
        std::optional<P3> best;
        double best_d = std::numeric_limits<double>::infinity();
        each([&](P3 p) {
            const double dd = distance(p, raw);
            if (dd <= tol && dd < best_d) {
                best = p;
                best_d = dd;
            }
        });
        return best;
    };
    if (const auto p = best_of([&](auto offer) {
            for (const Solid& s : solids) {
                for (const P3& v : s.vertices()) offer(v);
            }
        })) {
        return {*p, SnapKind::vertex};
    }
    const auto each_edge = [&](auto&& f) {
        for (const Solid& s : solids) {
            const auto v = s.vertices();
            for (const auto& [i, j] : s.edges()) f(v[static_cast<std::size_t>(i)], v[static_cast<std::size_t>(j)]);
        }
    };
    if (const auto p = best_of([&](auto offer) {
            each_edge([&](P3 a, P3 b) {
                offer({round_to((static_cast<double>(a.x) + static_cast<double>(b.x)) / 2), round_to((static_cast<double>(a.y) + static_cast<double>(b.y)) / 2),
                       round_to((static_cast<double>(a.z) + static_cast<double>(b.z)) / 2)});
            });
        })) {
        return {*p, SnapKind::midpoint};
    }
    if (const auto p = best_of([&](auto offer) {
            each_edge([&](P3 a, P3 b) {
                double t = 0.0, len2 = 0.0;
                for (int k = 0; k < 3; ++k) {
                    const double d = static_cast<double>(b[k] - a[k]);
                    t += d * static_cast<double>(raw[k] - a[k]);
                    len2 += d * d;
                }
                if (len2 == 0.0) return;
                t = std::clamp(t / len2, 0.0, 1.0);
                P3 q;
                for (int k = 0; k < 3; ++k) q[k] = round_to(static_cast<double>(a[k]) + t * static_cast<double>(b[k] - a[k]));
                // Along the edge on the grid too (Unit 26m slice 4): a pixel's place leaves slivers of a quantum.
                for (int k = 0; k < 3; ++k) {
                    if (grid > 0 && a[k] != b[k]) q[k] = std::clamp(on_grid(q[k]), std::min(a[k], b[k]), std::max(a[k], b[k]));
                }
                offer(q);
            });
        })) {
        return {*p, SnapKind::edge};
    }
    // Faces: the point onto a face's plane where its projection lies in the face; not the work plane's own plane.
    if (const auto p = best_of([&](auto offer) {
            const auto onto = [&](int a, Coord c, auto&& in_face) {
                if (a == plane.normal && c == plane.offset) return;
                P3 q = raw;
                q[a] = c;
                // In the face's plane on the grid too (slice 4), where that stays in the face.
                P3 g = q;
                for (int k = 0; k < 3; ++k) {
                    if (grid > 0 && k != a) g[k] = on_grid(g[k]);
                }
                if (in_face(g)) offer(g);
                else if (in_face(q)) offer(q);
            };
            for (const Solid& s : solids) {
                if (s.kind == SolidKind::box) {
                    for (int a = 0; a < 3; ++a) {
                        const auto in_face = [&](P3 q) {
                            const int u = (a + 1) % 3, v = (a + 2) % 3;
                            return q[u] >= s.box.lo[u] && q[u] <= s.box.hi[u] && q[v] >= s.box.lo[v] && q[v] <= s.box.hi[v];
                        };
                        onto(a, s.box.lo[a], in_face);
                        onto(a, s.box.hi[a], in_face);
                    }
                    continue;
                }
                const auto in_cap = [&](P3 q) { return inside(s.polygon, uv_of(s.axis, q)); };
                onto(s.axis, s.from, in_cap);
                onto(s.axis, s.to, in_cap);
                const std::size_t n = s.polygon.size();
                for (std::size_t k = 0; k < n; ++k) {
                    const Pt a = s.polygon[k], b = s.polygon[(k + 1) % n];
                    if (a.x != b.x && a.y != b.y) continue;
                    const bool along_u = a.y == b.y;  // the side's plane is v = a.y, else u = a.x
                    const int normal = along_u ? (s.axis + 2) % 3 : (s.axis + 1) % 3, along = along_u ? (s.axis + 1) % 3 : (s.axis + 2) % 3;
                    const Coord lo = along_u ? std::min(a.x, b.x) : std::min(a.y, b.y), hi = along_u ? std::max(a.x, b.x) : std::max(a.y, b.y);
                    onto(normal, along_u ? a.y : a.x, [&](P3 q) { return q[along] >= lo && q[along] <= hi && q[s.axis] >= s.from && q[s.axis] <= s.to; });
                }
            }
        })) {
        return {*p, SnapKind::face};
    }
    if (grid > 0) {
        const Pt uv = uv_of(plane.normal, raw);
        const auto g = [&](Coord c) { return round_to(static_cast<double>(c) / static_cast<double>(grid)) * grid; };
        return {from_uv(plane.normal, {g(uv.x), g(uv.y)}, plane.offset), SnapKind::grid};
    }
    return {raw, SnapKind::none};
}


// --- validation (slice 3, SD5) -------------------------------------------------------------------------------------

namespace {

// The points "(x, y, z)" in cm a message names, in quanta.
std::vector<P3> named_points(std::string_view m) {
    std::vector<P3> out;
    for (std::size_t open = m.find('('); open != std::string_view::npos; open = m.find('(', open + 1)) {
        const std::size_t close = m.find(')', open);
        if (close == std::string_view::npos) break;
        std::vector<double> v;
        std::string_view inner = m.substr(open + 1, close - open - 1);
        bool ok = true;
        while (ok && !inner.empty()) {
            const std::size_t comma = inner.find(',');
            const std::string t(trim(inner.substr(0, comma)));
            char* end = nullptr;
            const double x = t.empty() ? 0.0 : std::strtod(t.c_str(), &end);
            ok = !t.empty() && end == t.c_str() + t.size() && std::isfinite(x);
            if (ok) v.push_back(x);
            inner = comma == std::string_view::npos ? std::string_view{} : inner.substr(comma + 1);
        }
        if (ok && v.size() == 3) {
            out.push_back({static_cast<Coord>(std::llround(v[0] / geo::quantum_cm)), static_cast<Coord>(std::llround(v[1] / geo::quantum_cm)),
                           static_cast<Coord>(std::llround(v[2] / geo::quantum_cm))});
        }
    }
    return out;
}

}  // namespace

Problem locate_problem(const Model& m, std::string_view text, int line, std::string message) {
    Problem out{line, std::move(message), {}};
    const auto d = deck::parse(text);
    if (!d || line <= 0) return out;
    const deck::Section* at = nullptr;
    int region_index = -1, contact_index = -1, ri = 0, ci = 0;
    for (const auto& s : d->sections) {
        if (s.line > line) break;
        at = &s;
        region_index = s.kind == "region" ? ri : -1;
        contact_index = s.kind == "contact" ? ci : -1;
        ri += s.kind == "region" ? 1 : 0;
        ci += s.kind == "contact" ? 1 : 0;
    }
    if (!at) return out;
    const auto add = [&](const Marker& mk) {
        if (std::ranges::find(out.markers, mk) == out.markers.end()) out.markers.push_back(mk);
    };
    const auto contact_of = [&](int origin) {
        for (std::size_t k = 0; k < m.contacts.size(); ++k) {
            if (m.contacts[k].origin == origin) return static_cast<int>(k);
        }
        return -1;
    };
    if (region_index >= 0) {
        for (std::size_t r = 0; r < m.regions.size(); ++r) {
            if (m.regions[r].origin != region_index) continue;
            int solid = -1, n = 0;
            if (m.mesh == MeshKind::layered) {
                for (const auto& e : at->entries) {
                    if (!is_solid_key(e.key)) continue;
                    if (e.line == line) solid = n;
                    ++n;
                }
            }
            if (solid >= static_cast<int>(m.regions[r].solids.size())) solid = -1;
            add({m.regions[r].solids.empty() ? Marker::Kind::domain : Marker::Kind::solid, {static_cast<int>(r), solid}, -1, {}});
        }
    } else if (contact_index >= 0) {
        if (const int k = contact_of(contact_index); k >= 0) add({Marker::Kind::contact, {}, k, {}});
    } else if (at->kind == "mesh") {
        const deck::Entry* patch = nullptr;
        for (const auto& e : at->entries) {
            if (e.line == line && e.key.starts_with("patch ")) patch = &e;
        }
        if (patch) {
            const std::string_view name = trim(std::string_view(patch->key).substr(6));
            int c = 0;
            for (const auto& s : d->sections) {
                if (s.kind != "contact") continue;
                const deck::Entry* p = find_entry(s, "patch");
                if (p && trim(p->value) == name) {
                    if (const int k = contact_of(c); k >= 0) add({Marker::Kind::contact, {}, k, {}});
                }
                ++c;
            }
        }
    } else {
        return out;  // outside the geometry
    }
    const auto points = named_points(out.message);
    for (const P3 p : points) {
        add({Marker::Kind::point, {}, -1, p});
        for (std::size_t r = 0; r < m.regions.size(); ++r) {
            for (std::size_t k = 0; k < m.regions[r].solids.size(); ++k) {
                const Box3 b = m.regions[r].solids[k].bounds();
                bool in = true;
                for (int a = 0; a < 3; ++a) in = in && p[a] >= b.lo[a] - 1 && p[a] <= b.hi[a] + 1;
                if (in) add({Marker::Kind::solid, {static_cast<int>(r), static_cast<int>(k)}, -1, {}});
            }
        }
    }
    if (at->kind == "mesh" && out.markers.empty()) add({Marker::Kind::domain, {}, -1, {}});
    return out;
}


int region_at(const Model& m, P3 p) {
    const auto holds = [&](const Region& reg) {
        return std::ranges::any_of(reg.solids, [&](const Solid& s) {
            if (s.kind == SolidKind::box) {
                return p.x > s.box.lo.x && p.x < s.box.hi.x && p.y > s.box.lo.y && p.y < s.box.hi.y && p.z > s.box.lo.z && p.z < s.box.hi.z;
            }
            return p[s.axis] > s.from && p[s.axis] < s.to && inside(s.polygon, uv_of(s.axis, p));
        });
    };
    for (const int r : std::views::reverse(paint_order(m))) {
        if (holds(m.regions[static_cast<std::size_t>(r)])) return r;
    }
    return m.background();
}

std::vector<int> paint_order(const Model& m) {
    std::vector<int> order;
    for (std::size_t r = m.regions.size(); r-- > 0;) {
        if (m.regions[r].fill && !m.regions[r].solids.empty()) order.push_back(static_cast<int>(r));
    }
    for (std::size_t r = 0; r < m.regions.size(); ++r) {
        if (!m.regions[r].fill && !m.regions[r].solids.empty()) order.push_back(static_cast<int>(r));
    }
    return order;
}

// --- parameters (Unit 26m MT8) -----------------------------------------------------------------------------------

Coord* slot(Solid& s, int k) {
    if (k < 0) return nullptr;
    if (s.kind == SolidKind::box) {
        Coord* const c[6] = {&s.box.lo.x, &s.box.hi.x, &s.box.lo.y, &s.box.hi.y, &s.box.lo.z, &s.box.hi.z};
        return k < 6 ? c[k] : nullptr;
    }
    if (k == 0) return &s.from;
    if (k == 1) return &s.to;
    const auto j = static_cast<std::size_t>((k - 2) / 2);
    if (j >= s.polygon.size()) return nullptr;
    return (k - 2) % 2 == 0 ? &s.polygon[j].x : &s.polygon[j].y;
}

int slot_count(const Solid& s) { return s.kind == SolidKind::box ? 6 : 2 + 2 * static_cast<int>(s.polygon.size()); }

namespace {

Coord* box_slot(Box3& b, int k) {
    Coord* const c[6] = {&b.lo.x, &b.hi.x, &b.lo.y, &b.hi.y, &b.lo.z, &b.hi.z};
    return k >= 0 && k < 6 ? c[k] : nullptr;
}

bool param_name(std::string_view n) {
    if (n.empty() || !(std::isalpha(static_cast<unsigned char>(n.front())) || n.front() == '_')) return false;
    return std::ranges::all_of(n, [](char c) { return std::isalnum(static_cast<unsigned char>(c)) || c == '_'; });
}

// A tied coordinate: where it is, and its reference.
struct Tie {
    Coord* at = nullptr;
    std::string* ref = nullptr;
};

std::vector<Tie> ties(Model& m) {
    std::vector<Tie> out;
    for (Region& r : m.regions) {
        for (Solid& s : r.solids) {
            for (std::size_t k = 0; k < s.refs.size(); ++k) {
                if (Coord* c = slot(s, static_cast<int>(k)); c && !s.refs[k].empty()) out.push_back({c, &s.refs[k]});
            }
        }
    }
    for (std::size_t k = 0; k < m.domain_refs.size(); ++k) {
        if (Coord* c = box_slot(m.domain, static_cast<int>(k)); c && !m.domain_refs[k].empty()) out.push_back({c, &m.domain_refs[k]});
    }
    return out;
}

// Every parameter resolved again from its expression, then every tied coordinate; the solids and domain checked.
std::expected<void, std::string> resolve_all(Model& m) {
    std::vector<deck::Entry> entries;
    for (std::size_t k = 0; k < m.parameters.size(); ++k) entries.push_back({m.parameters[k].name, m.parameters[k].expr, static_cast<int>(k) + 1});
    const auto r = deck::resolve_parameters(entries, "[parameters]");
    if (!r) return std::unexpected(r.error().message);
    for (Param& p : m.parameters) {
        const auto it = std::ranges::find(*r, p.name, &deck::Parameter::name);
        p.text = it->text;
        p.length = it->length;
        p.value = 0;
        p.number = 0;
        if (p.length) {
            const auto q = parse_length(it->text);
            if (!q) return std::unexpected("parameter '" + p.name + "' is not a length the modeller can draw");
            p.value = *q;
        } else {
            p.number = it->value;
        }
    }
    const Params dp = deck_params(m);
    for (const Tie& t : ties(m)) {
        const auto v = resolved(*t.ref, dp);
        if (!v) return std::unexpected("'" + *t.ref + "' does not give a length");
        *t.at = *v;
    }
    for (const Region& reg : m.regions) {
        for (const Solid& s : reg.solids) {
            if (!solid_ok(s)) return std::unexpected(std::string("a tied solid would not be valid (a face past its opposite)"));
        }
    }
    if (!box_ok(m.domain)) return std::unexpected(std::string("the tied domain would not be valid"));
    return {};
}

Param* find_param(Model& m, std::string_view name) {
    const auto it = std::ranges::find(m.parameters, name, &Param::name);
    return it == m.parameters.end() ? nullptr : &*it;
}

// A lone reference `$p`: p; else none.
std::optional<std::string_view> lone(std::string_view ref) {
    ref = trim(ref);
    if (ref.size() < 2 || ref.front() != '$' || !param_name(ref.substr(1))) return std::nullopt;
    return ref.substr(1);
}

// A sum's last parameter and its sign (`... + $q`, `... - $q`; a leading `$q` is +).
std::optional<std::pair<std::string_view, int>> last_term(std::string_view ref) {
    const std::size_t at = ref.rfind('$');
    if (at == std::string_view::npos) return std::nullopt;
    std::size_t end = at + 1;
    while (end < ref.size() && (std::isalnum(static_cast<unsigned char>(ref[end])) || ref[end] == '_')) ++end;
    std::size_t k = at;
    while (k > 0 && (ref[k - 1] == ' ' || ref[k - 1] == '\t')) --k;
    const int sign = k > 0 && ref[k - 1] == '-' ? -1 : 1;
    return std::pair{ref.substr(at + 1, end - at - 1), sign};
}

// A parameter a drag may set: a length written as a value.
bool settable(const Param* p) { return p && p->length && parse_length(p->expr).has_value(); }

// A length for a parameter's line: always with a unit (a bare 0 would be a number, which no sum takes).
std::string length_text(Coord v) {
    std::string s = format_length(v);
    if (!s.ends_with("nm") && !s.ends_with("um") && !s.ends_with("cm")) s += "nm";
    return s;
}

}  // namespace

std::expected<void, std::string> set_parameter(Model& m, std::string_view name, std::string expr) {
    if (!param_name(name)) return std::unexpected("'" + std::string(name) + "' is not a parameter name (a letter or _, then letters, digits, _)");
    Model t = m;
    if (Param* p = find_param(t, name)) {
        p->expr = std::move(expr);
    } else {
        t.parameters.push_back({std::string(name), std::move(expr), {}, false, 0, 0});
    }
    if (auto r = resolve_all(t); !r) return r;
    m = std::move(t);
    return {};
}

void tie_edits(const Model& before, Model& after) {
    if (before.parameters != after.parameters) return;  // the parameters themselves edited: nothing to follow
    Model t = after;
    struct Moved {
        Coord* at;
        std::string* ref;
        Coord want;
    };
    std::vector<Moved> moved;
    const auto collect = [&](const std::vector<std::string>& was_refs, std::vector<std::string>& refs, auto&& was_at, auto&& at, int n) {
        if (was_refs != refs) return;
        bool edited = false;  // any coordinate of it moved: its tied ones are all held where the edit put them
        for (int k = 0; k < n; ++k) edited = edited || *at(k) != *was_at(k);
        if (!edited) return;
        for (const bool changed : {true, false}) {  // the moved ones first: theirs win a parameter both set
            for (int k = 0; k < n && k < static_cast<int>(refs.size()); ++k) {
                if (refs[static_cast<std::size_t>(k)].empty()) continue;
                Coord* now = at(k);
                if ((*now != *was_at(k)) == changed) moved.push_back({now, &refs[static_cast<std::size_t>(k)], *now});
            }
        }
    };
    for (std::size_t i = 0; i < t.regions.size() && i < before.regions.size(); ++i) {
        for (std::size_t j = 0; j < t.regions[i].solids.size() && j < before.regions[i].solids.size(); ++j) {
            Solid& a = t.regions[i].solids[j];
            Solid b = before.regions[i].solids[j];
            if (slot_count(a) != slot_count(b) || a.kind != b.kind) continue;
            collect(b.refs, a.refs, [&](int k) { return slot(b, k); }, [&](int k) { return slot(a, k); }, slot_count(a));
        }
    }
    Box3 was_domain = before.domain;
    collect(before.domain_refs, t.domain_refs, [&](int k) { return box_slot(was_domain, k); }, [&](int k) { return box_slot(t.domain, k); }, 6);
    if (moved.empty()) return;

    std::map<std::string, Coord, std::less<>> set;  // the parameters set, by name
    const auto assign = [&](std::string_view name, Coord v) {
        Param* p = find_param(t, name);
        if (!settable(p)) return false;
        if (const auto it = set.find(name); it != set.end()) return it->second == v;
        set.emplace(std::string(name), v);
        p->expr = length_text(v);
        return true;
    };
    // Lone references first: each sets its parameter to where its coordinate went.
    for (Moved& c : moved) {
        if (const auto n = lone(*c.ref); n && !assign(*n, c.want)) c.ref->clear();
    }
    bool ok = static_cast<bool>(resolve_all(t));
    // Then sums still off: the last parameter takes the rest.
    if (ok) {
        const Params dp = deck_params(t);
        for (Moved& c : moved) {
            if (c.ref->empty() || lone(*c.ref)) continue;
            const auto got = resolved(*c.ref, dp);
            if (got == c.want) continue;
            const auto q = last_term(*c.ref);
            const Param* p = q ? find_param(t, q->first) : nullptr;
            if (!got || !settable(p) || !assign(q->first, p->value + q->second * (c.want - *got))) c.ref->clear();
        }
        ok = static_cast<bool>(resolve_all(t));
    }
    if (!ok) {
        // Nothing resolves as set: every moved tie broken, the parameters as they were, the edit's values kept.
        t = after;
        moved.clear();
        for (std::size_t i = 0; i < t.regions.size() && i < before.regions.size(); ++i) {
            for (std::size_t j = 0; j < t.regions[i].solids.size() && j < before.regions[i].solids.size(); ++j) {
                Solid& a = t.regions[i].solids[j];
                Solid b = before.regions[i].solids[j];
                if (slot_count(a) != slot_count(b) || a.kind != b.kind) continue;
                collect(b.refs, a.refs, [&](int k) { return slot(b, k); }, [&](int k) { return slot(a, k); }, slot_count(a));
            }
        }
        collect(before.domain_refs, t.domain_refs, [&](int k) { return box_slot(was_domain, k); }, [&](int k) { return box_slot(t.domain, k); }, 6);
        for (Moved& c : moved) c.ref->clear();
        after = std::move(t);
        return;
    }
    // A moved coordinate its tie does not give: the tie broken, the edit's value kept.
    for (Moved& c : moved) {
        if (*c.at != c.want) {
            c.ref->clear();
            *c.at = c.want;
        }
    }
    after = std::move(t);
}

bool name_extent(Model& m, SolidRef s, int near, int far, std::string_view name) {
    if (!m.editable() || !valid(m, s) || !param_name(name) || near == far) return false;
    Model t = m;
    Solid& solid = t.regions[static_cast<std::size_t>(s.region)].solids[static_cast<std::size_t>(s.solid)];
    Coord* a = slot(solid, near);
    Coord* b = slot(solid, far);
    if (!a || !b) return false;
    solid.refs.resize(static_cast<std::size_t>(slot_count(solid)));
    const std::string near_text = solid.refs[static_cast<std::size_t>(near)].empty() ? format_length(*a) : solid.refs[static_cast<std::size_t>(near)];
    const Param* p = find_param(t, name);
    if (p && !p->length) return false;
    const bool down = far < near;  // the low face (a box's lo, a prism's from) lies below its opposite: near - $name
    if (!p) t.parameters.push_back({std::string(name), length_text(down ? *a - *b : *b - *a), {}, false, 0, 0});
    solid.refs[static_cast<std::size_t>(far)] = near_text + (down ? " - $" : " + $") + std::string(name);
    if (!resolve_all(t)) return false;
    m = std::move(t);
    return true;
}

// --- doping and mesh tools (Unit 26m slice 4: MT9, MT10) ----------------------------------------------------------

namespace {

double cm(Coord q) { return static_cast<double>(q) / 1e8; }  // correctly rounded (q * 1e-8 is not: 7000 gives 7.000000000000001e-05)

// `line` added as the last entry of the deck's `kind` section; a new section (before the first contact or analysis
// section, else at the end) when it has none. Every other line kept byte for byte.
std::expected<std::string, std::string> append_entry(std::string_view text, std::string_view kind, const std::string& line) {
    const auto d = deck::parse(text);
    if (!d) return std::unexpected(d.error().message);
    Lines l = split(text);
    const deck::Section* section = nullptr;
    for (const auto& s : d->sections) {
        if (s.kind == kind) section = &s;
    }
    if (section) {
        const int at = section_last(l, *section) + 1;
        l.body.insert(l.body.begin() + at, line);
        l.eol.insert(l.eol.begin() + at, l.nl);
        if (at > 0 && l.eol[static_cast<std::size_t>(at - 1)].empty()) l.eol[static_cast<std::size_t>(at - 1)] = l.nl;
    } else {
        static const std::set<std::string_view> later{"contact", "thermal_contact", "equilibrium", "sweep", "transient",
                                                      "small_signal", "trace", "output", "adapt"};
        int at = static_cast<int>(l.body.size());
        for (const auto& s : d->sections) {
            if (later.contains(s.kind)) {
                at = s.line - 1;
                break;
            }
        }
        const std::vector<std::string> block{"[" + std::string(kind) + "]", line, ""};
        if (at > 0 && l.eol[static_cast<std::size_t>(at - 1)].empty()) l.eol[static_cast<std::size_t>(at - 1)] = l.nl;
        for (std::size_t k = 0; k < block.size(); ++k) {
            l.body.insert(l.body.begin() + at + static_cast<int>(k), block[k]);
            l.eol.insert(l.eol.begin() + at + static_cast<int>(k), l.nl);
        }
    }
    std::string out;
    for (std::size_t k = 0; k < l.body.size(); ++k) out += l.body[k] + l.eol[k];
    return out;
}

// erfc⁻¹(y) for 0 < y < 2, by Newton's method on erfc from a rational start (to the last bits).
double erfc_inverse(double y) {
    double x = y < 1.0 ? std::sqrt(-std::log(y / 2.0)) : -std::sqrt(-std::log(1.0 - y / 2.0));
    for (int k = 0; k < 60; ++k) {
        const double f = std::erfc(x) - y;
        const double step = f / (-2.0 / std::sqrt(std::numbers::pi) * std::exp(-x * x));
        x -= step;
        if (std::abs(step) <= 1e-16 * std::max(1.0, std::abs(x))) break;
    }
    return x;
}

}  // namespace

std::vector<device::DopingProfile> profiles(std::string_view text) {
    std::vector<device::DopingProfile> out;
    const auto d = deck::parse(text);
    if (!d) return out;
    const auto params = deck::parameters(*d);
    for (const auto& s : d->sections) {
        if (s.kind != "doping") continue;
        for (const auto& e : s.entries) {
            if (e.key != "profile") continue;
            const auto v = params ? deck::substitute(e.value, *params) : std::expected<std::string, std::string>(e.value);
            if (!v) continue;
            if (auto p = deck::parse_profile(*v)) out.push_back(*p);
        }
    }
    return out;
}

std::expected<std::string, std::string> add_profile(std::string_view text, const device::DopingProfile& p) {
    return append_entry(text, "doping", "profile = " + deck::profile_text(p));
}

device::DopingProfile region_profile(const Model& m, SolidRef s, device::Dopant dopant, double peak) {
    device::DopingProfile p{.dopant = dopant, .peak_cm3 = peak, .shape = device::DepthShape::constant};
    const Box3 b = m.regions[static_cast<std::size_t>(s.region)].solids[static_cast<std::size_t>(s.solid)].bounds();
    // The depth axis (unbounded: the deck's windows are on the other two) the one the solid spans most of the domain
    // along, z on a tie: exact when it spans the domain there.
    double best = -1.0;
    for (const int a : {2, 1, 0}) {
        const double span = static_cast<double>(b.hi[a] - b.lo[a]) / static_cast<double>(std::max<Coord>(1, m.domain.hi[a] - m.domain.lo[a]));
        if (span > best + 1e-12) {
            best = span;
            p.axis = a;
        }
    }
    for (int a = 0; a < 3; ++a) {
        if (a == p.axis) continue;
        p.lo[static_cast<std::size_t>(a)] = cm(b.lo[a]);
        p.hi[static_cast<std::size_t>(a)] = cm(b.hi[a]);
    }
    return p;
}

device::DopingProfile implant_profile(Normal face, Coord at_face, const Box3& window, Coord depth, Coord lateral, device::Dopant dopant,
                                      double peak) {
    device::DopingProfile p{.dopant = dopant, .peak_cm3 = peak, .shape = device::DepthShape::gaussian, .axis = face.axis};
    p.center = cm(at_face - face.sign * depth);  // into the solid: against the face's outward normal
    p.width = cm(depth) / 2.0;
    p.lateral = cm(lateral);
    for (int a = 0; a < 3; ++a) {
        if (a == face.axis) continue;
        p.lo[static_cast<std::size_t>(a)] = cm(window.lo[a]);
        p.hi[static_cast<std::size_t>(a)] = cm(window.hi[a]);
    }
    return p;
}

std::expected<device::DopingProfile, std::string> junction_profile(const std::vector<device::DopingProfile>& existing, Normal face,
                                                                   Coord at_face, const Box3& window, Coord depth,
                                                                   device::Dopant dopant, double peak) {
    if (face.sign > 0) {
        return std::unexpected(std::string("a junction diffuses from a face looking down an axis (-x, -y or -z): the erfc profile falls along its axis"));
    }
    if (depth <= 0) return std::unexpected(std::string("the junction lies inside the solid: drag into it"));
    device::DopingProfile p{.dopant = dopant, .peak_cm3 = peak, .shape = device::DepthShape::erfc, .axis = face.axis};
    for (int a = 0; a < 3; ++a) {
        if (a == face.axis) continue;
        p.lo[static_cast<std::size_t>(a)] = cm(window.lo[a]);
        p.hi[static_cast<std::size_t>(a)] = cm(window.hi[a]);
    }
    // The junction's point: the window's centre at the dragged depth; the background there the other dopant's net.
    mesh::Point x{};
    for (int a = 0; a < 3; ++a) x[static_cast<std::size_t>(a)] = 0.5 * (cm(window.lo[a]) + cm(window.hi[a]));
    x[static_cast<std::size_t>(face.axis)] = cm(at_face + depth);
    double background = 0.0;
    for (const auto& q : existing) background += (q.dopant == dopant ? -1.0 : 1.0) * device::evaluate(q, x);
    if (!(background > 0.0)) return std::unexpected(std::string("no opposite doping at that depth for a junction to meet"));
    if (!(background < peak)) return std::unexpected(std::string("the opposite doping there is above the level: raise the level"));
    // peak/2 erfc((xj - c)/w) = background, with w half the depth: c = xj - w erfc⁻¹(2 background / peak).
    p.width = cm(depth) / 2.0;
    p.center = x[static_cast<std::size_t>(face.axis)] - p.width * erfc_inverse(2.0 * background / peak);
    return p;
}


Box3 profile_extent(const Box3& domain, const device::DopingProfile& p) {
    const auto q = [](double c) { return static_cast<Coord>(std::llround(c * 1e8)); };
    Box3 b = domain;
    for (int a = 0; a < 3; ++a) {
        const auto i = static_cast<std::size_t>(a);
        double lo = -std::numeric_limits<double>::infinity(), hi = std::numeric_limits<double>::infinity();
        if (a != p.axis || p.shape == device::DepthShape::constant) {  // Unit 67: a constant is bounded on every axis
            lo = p.lo[i] - 2.0 * p.lateral;
            hi = p.hi[i] + 2.0 * p.lateral;
        } else if (p.shape == device::DepthShape::gaussian) {
            lo = p.center - 2.0 * p.width;
            hi = p.center + 2.0 * p.width;
        } else if (p.shape == device::DepthShape::erfc) {
            hi = p.center + 2.0 * p.width;
        }
        if (std::isfinite(lo)) b.lo[a] = std::clamp(q(lo), domain.lo[a], domain.hi[a]);
        if (std::isfinite(hi)) b.hi[a] = std::clamp(q(hi), domain.lo[a], domain.hi[a]);
        b.lo[a] = std::min(b.lo[a], b.hi[a]);
    }
    return b;
}
std::vector<SizeBox3> size_boxes(std::string_view text) {
    std::vector<SizeBox3> out;
    const auto d = deck::parse(text);
    if (!d) return out;
    for (const auto& s : d->sections) {
        if (s.kind != "mesh") continue;
        for (const auto& e : s.entries) {
            if (e.key != "box") continue;
            std::vector<std::string> refs;
            const auto c = lengths(e.value, 7, {}, refs);
            if (!c) continue;
            out.push_back({normalized({(*c)[0], (*c)[1], (*c)[2]}, {(*c)[3], (*c)[4], (*c)[5]}), (*c)[6]});
        }
    }
    return out;
}

std::expected<std::string, std::string> add_size_box(std::string_view text, const SizeBox3& b) {
    if (b.size <= 0 || !(b.box.lo.x < b.box.hi.x && b.box.lo.y < b.box.hi.y && b.box.lo.z < b.box.hi.z)) {
        return std::unexpected(std::string("a refine box needs a volume and a size"));
    }
    return append_entry(text, "mesh", "box = " + patch_text(b.box) + "  " + format_length(b.size));
}

}  // namespace NiTCAD::app::solid
