// The 3D solid modeller's model (solid.hpp; DECISIONS.md Unit 26l, SD3).
#include "NiTCAD/app/solid.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>

#include "deck_text.hpp"
#include "NiTCAD/deck/syntax.hpp"

namespace NiTCAD::app::solid {

namespace {

using namespace deck_text;
using geo::format_length;
using geo::parse_length;

std::optional<std::vector<Coord>> lengths(std::string_view value, std::size_t n) {
    const auto w = words(value);
    if (w.size() != n) return std::nullopt;
    std::vector<Coord> out;
    for (const auto t : w) {
        const auto v = parse_length(t);
        if (!v) return std::nullopt;
        out.push_back(*v);
    }
    return out;
}

// A region box or the domain: x0 x1 y0 y1 z0 z1.
std::optional<Box3> parse_box3(std::string_view value) {
    const auto c = lengths(value, 6);
    if (!c || (*c)[0] > (*c)[1] || (*c)[2] > (*c)[3] || (*c)[4] > (*c)[5]) return std::nullopt;
    return Box3{{(*c)[0], (*c)[2], (*c)[4]}, {(*c)[1], (*c)[3], (*c)[5]}};
}

// A patch box: x0 y0 z0 x1 y1 z1.
std::optional<Box3> parse_patch(std::string_view value) {
    const auto c = lengths(value, 6);
    if (!c || (*c)[0] > (*c)[3] || (*c)[1] > (*c)[4] || (*c)[2] > (*c)[5]) return std::nullopt;
    return Box3{{(*c)[0], (*c)[1], (*c)[2]}, {(*c)[3], (*c)[4], (*c)[5]}};
}

std::optional<std::vector<Pt>> parse_points(std::string_view value) {
    const auto w = words(value);
    if (w.size() < 6 || w.size() % 2 != 0) return std::nullopt;
    std::vector<Pt> out;
    for (std::size_t k = 0; k < w.size(); k += 2) {
        const auto u = parse_length(w[k]), v = parse_length(w[k + 1]);
        if (!u || !v) return std::nullopt;
        out.push_back({*u, *v});
    }
    return out;
}

std::optional<int> axis_of(std::string_view a) {
    if (a == "x") return 0;
    if (a == "y") return 1;
    if (a == "z") return 2;
    return std::nullopt;
}

std::optional<Solid> parse_solid(const deck::Entry& e) {
    if (e.key == "box") {
        const auto b = parse_box3(e.value);
        if (!b) return std::nullopt;
        return Solid{SolidKind::box, *b, 2, 0, 0, {}};
    }
    if (e.key == "prism") {
        // <axis> <from> <to>, u v, u v, ...
        const std::string_view text = e.value;
        const std::size_t comma = text.find(',');
        if (comma == std::string_view::npos) return std::nullopt;
        const auto head = words(text.substr(0, comma));
        if (head.size() != 3) return std::nullopt;
        const auto axis = axis_of(head[0]);
        const auto from = parse_length(head[1]), to = parse_length(head[2]);
        auto p = parse_points(text.substr(comma + 1));
        if (!axis || !from || !to || !p) return std::nullopt;
        Solid s;
        s.kind = SolidKind::prism;
        s.axis = *axis;
        s.from = *from;
        s.to = *to;
        s.polygon = std::move(*p);
        return s;
    }
    return std::nullopt;
}

bool is_solid_key(std::string_view k) { return k == "box" || k == "prism"; }

constexpr const char* axis_name[] = {"x", "y", "z"};

std::string box_text(const Box3& b) {  // x0 x1 y0 y1 z0 z1
    return format_length(b.lo.x) + " " + format_length(b.hi.x) + " " + format_length(b.lo.y) + " " + format_length(b.hi.y) + " " +
           format_length(b.lo.z) + " " + format_length(b.hi.z);
}

std::string patch_text(const Box3& b) {  // x0 y0 z0  x1 y1 z1
    return format_length(b.lo.x) + " " + format_length(b.lo.y) + " " + format_length(b.lo.z) + "  " + format_length(b.hi.x) + " " +
           format_length(b.hi.y) + " " + format_length(b.hi.z);
}

std::string solid_line(const Solid& s) {
    if (s.kind == SolidKind::box) return "box = " + box_text(s.box);
    std::string out = std::string("prism = ") + axis_name[s.axis] + " " + format_length(s.from) + " " + format_length(s.to);
    for (const Pt& p : s.polygon) out += ", " + format_length(p.x) + " " + format_length(p.y);
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
            const auto b = parse_patch(e.value);
            if (w.size() != 2 || !b) return refuse(at_line(e, "a patch box the modeller cannot read"));
            patches[std::string(w[1])] = *b;
        } else if (m.mesh == MeshKind::layered && e.key == "domain") {
            const auto b = parse_box3(e.value);
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
            if (m.mesh == MeshKind::layered && is_solid_key(e.key)) {
                auto solid = parse_solid(e);
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
            for (const Solid& s : r.solids) out.push_back(solid_line(s));
            return out;
        }
        const deck::Section& s = *regions[static_cast<std::size_t>(r.origin)];
        SectionEdit e;
        if (r.name != s.name) e.header = "[region " + r.name + "]";
        std::size_t solid = 0;
        int last_solid = -1;
        const int last_entry = s.entries.empty() ? -1 : static_cast<int>(s.entries.size()) - 1;
        bool material = false;
        for (std::size_t k = 0; k < s.entries.size(); ++k) {
            const auto& en = s.entries[k];
            if (en.key == "material") {
                material = true;
                if (std::string(trim(en.value)) != r.material) {
                    e.replace[k] = r.material.empty() ? std::nullopt : std::optional<std::string>("material = " + r.material);
                }
            } else if (solid_entry(en.key)) {
                last_solid = static_cast<int>(k);
                if (solid >= r.solids.size()) {
                    e.replace[k] = std::nullopt;
                } else if (en.key == "extrude" || parse_solid(en) != r.solids[solid]) {
                    e.replace[k] = solid_line(r.solids[solid]);
                }
                ++solid;
            }
        }
        if (!material && !r.material.empty()) e.after.emplace(-1, "material = " + r.material);
        for (; solid < r.solids.size(); ++solid) e.after.emplace(last_solid >= 0 ? last_solid : last_entry, solid_line(r.solids[solid]));
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
        const std::string domain_line = "domain = " + box_text(m.domain);
        if (domain_entry >= 0) {
            if (parse_box3(mesh->entries[static_cast<std::size_t>(domain_entry)].value) != m.domain) e.replace[static_cast<std::size_t>(domain_entry)] = domain_line;
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
                if (!used.contains(pname) && parse_patch(mesh->entries[it->second].value) != c.box) e.replace[it->second] = line;
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
                if (in_face(q)) offer(q);
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
    for (std::size_t r = m.regions.size(); r-- > 0;) {
        for (const Solid& s : m.regions[r].solids) {
            if (s.kind == SolidKind::box) {
                if (p.x > s.box.lo.x && p.x < s.box.hi.x && p.y > s.box.lo.y && p.y < s.box.hi.y && p.z > s.box.lo.z && p.z < s.box.hi.z) {
                    return static_cast<int>(r);
                }
            } else if (p[s.axis] > s.from && p[s.axis] < s.to && inside(s.polygon, uv_of(s.axis, p))) {
                return static_cast<int>(r);
            }
        }
    }
    return m.background();
}

}  // namespace NiTCAD::app::solid
