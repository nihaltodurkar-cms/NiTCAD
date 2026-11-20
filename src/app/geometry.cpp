// The geometry editor's model (geometry.hpp; DECISIONS.md Unit 26k, GD3-GD5).
#include "NiTCAD/app/geometry.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <format>
#include <map>
#include <set>
#include <utility>

#include "NiTCAD/deck/syntax.hpp"

namespace NiTCAD::app::geo {

namespace {

bool blank(char c) { return c == ' ' || c == '\t' || c == '\r'; }

std::string_view trim(std::string_view s) {
    while (!s.empty() && blank(s.front())) s.remove_prefix(1);
    while (!s.empty() && blank(s.back())) s.remove_suffix(1);
    return s;
}

// Tokens separated by blanks and commas; `groups` the comma-separated groups' tokens.
std::vector<std::string_view> words(std::string_view s) {
    std::vector<std::string_view> out;
    std::size_t i = 0;
    while (i < s.size()) {
        while (i < s.size() && (blank(s[i]) || s[i] == ',')) ++i;
        std::size_t j = i;
        while (j < s.size() && !blank(s[j]) && s[j] != ',') ++j;
        if (j > i) out.push_back(s.substr(i, j - i));
        i = j;
    }
    return out;
}

std::optional<Rect> parse_box(std::string_view value) {
    const auto w = words(value);
    if (w.size() != 4) return std::nullopt;
    Coord c[4];
    for (int k = 0; k < 4; ++k) {
        const auto v = parse_length(w[static_cast<std::size_t>(k)]);
        if (!v) return std::nullopt;
        c[k] = *v;
    }
    if (c[0] > c[1] || c[2] > c[3]) return std::nullopt;
    return Rect{{c[0], c[2]}, {c[1], c[3]}};
}

std::optional<std::vector<Pt>> parse_points(std::string_view value) {
    const auto w = words(value);
    if (w.size() < 2 || w.size() % 2 != 0) return std::nullopt;
    std::vector<Pt> out;
    for (std::size_t k = 0; k < w.size(); k += 2) {
        const auto x = parse_length(w[k]), y = parse_length(w[k + 1]);
        if (!x || !y) return std::nullopt;
        out.push_back({*x, *y});
    }
    return out;
}

std::optional<Shape> parse_shape(const deck::Entry& e) {
    if (e.key == "box") {
        const auto r = parse_box(e.value);
        if (!r) return std::nullopt;
        return Shape{ShapeKind::box, {r->lo, r->hi}};
    }
    if (e.key == "polygon") {
        auto p = parse_points(e.value);
        if (!p || p->size() < 3) return std::nullopt;
        return Shape{ShapeKind::polygon, std::move(*p)};
    }
    return std::nullopt;
}

bool is_shape_key(std::string_view k) { return k == "box" || k == "polygon"; }

std::string shape_line(const Shape& s) {
    if (s.kind == ShapeKind::box) {
        return "box = " + format_length(s.points[0].x) + " " + format_length(s.points[1].x) + " " +
               format_length(s.points[0].y) + " " + format_length(s.points[1].y);
    }
    std::string out = "polygon = ";
    for (std::size_t k = 0; k < s.points.size(); ++k) {
        if (k > 0) out += ", ";
        out += format_length(s.points[k].x) + " " + format_length(s.points[k].y);
    }
    return out;
}

std::string segment_text(Pt a, Pt b) {
    return format_length(a.x) + " " + format_length(a.y) + ", " + format_length(b.x) + " " + format_length(b.y);
}

const deck::Entry* find_entry(const deck::Section& s, std::string_view key) {
    for (const auto& e : s.entries) {
        if (e.key == key) return &e;
    }
    return nullptr;
}

std::optional<bool> on_off(std::string_view v) {
    v = trim(v);
    if (v == "on" || v == "yes" || v == "true") return true;
    if (v == "off" || v == "no" || v == "false") return false;
    return std::nullopt;
}

// A tensor deck's face: the side of the domain.
std::optional<std::pair<Pt, Pt>> face_segment(const Rect& d, std::string_view face) {
    if (face == "x_min") return std::pair{Pt{d.lo.x, d.lo.y}, Pt{d.lo.x, d.hi.y}};
    if (face == "x_max") return std::pair{Pt{d.hi.x, d.lo.y}, Pt{d.hi.x, d.hi.y}};
    if (face == "y_min") return std::pair{Pt{d.lo.x, d.lo.y}, Pt{d.hi.x, d.lo.y}};
    if (face == "y_max") return std::pair{Pt{d.lo.x, d.hi.y}, Pt{d.hi.x, d.hi.y}};
    return std::nullopt;
}

std::string face_of(const Rect& d, Pt a, Pt b) {
    if (a.x == b.x && a.x == d.lo.x) return "x_min";
    if (a.x == b.x && a.x == d.hi.x) return "x_max";
    if (a.y == b.y && a.y == d.lo.y) return "y_min";
    return "y_max";
}

// --- the text as lines ----------------------------------------------------------------------------------------

struct Lines {
    std::vector<std::string> body, eol;
    std::string nl = "\n";
};

Lines split(std::string_view t) {
    Lines l;
    if (t.find("\r\n") != std::string_view::npos) l.nl = "\r\n";
    std::size_t i = 0;
    while (i < t.size()) {
        const std::size_t j = t.find('\n', i);
        if (j == std::string_view::npos) {
            l.body.emplace_back(t.substr(i));
            l.eol.emplace_back();
            break;
        }
        const bool cr = j > i && t[j - 1] == '\r';
        l.body.emplace_back(t.substr(i, j - i - (cr ? 1 : 0)));
        l.eol.emplace_back(cr ? "\r\n" : "\n");
        i = j + 1;
    }
    return l;
}

bool continues(std::string_view s) {
    const std::size_t hash = s.find('#');
    if (hash != std::string_view::npos) s = s.substr(0, hash);
    s = trim(s);
    return !s.empty() && s.back() == '\\';
}

// 0-based line spans.
int entry_last(const Lines& l, int first) {
    int k = first;
    while (k + 1 < static_cast<int>(l.body.size()) && continues(l.body[static_cast<std::size_t>(k)])) ++k;
    return k;
}

int section_last(const Lines& l, const deck::Section& s) {
    int last = s.line - 1;
    for (const auto& e : s.entries) last = std::max(last, entry_last(l, e.line - 1));
    return last;
}

struct Change {
    int first = 0, count = 0;  // the lines replaced (count 0: an insertion before `first`)
    std::vector<std::string> lines;
};

// A section's block rewritten: `edit` sees its lines and the entries' spans relative to the block.
struct Block {
    std::vector<std::string> lines;
};

}  // namespace

Rect normalized(Pt a, Pt b) noexcept { return {{std::min(a.x, b.x), std::min(a.y, b.y)}, {std::max(a.x, b.x), std::max(a.y, b.y)}}; }

Rect Shape::bounds() const noexcept {
    Rect r{points.front(), points.front()};
    for (const Pt& p : points) {
        r.lo = {std::min(r.lo.x, p.x), std::min(r.lo.y, p.y)};
        r.hi = {std::max(r.hi.x, p.x), std::max(r.hi.y, p.y)};
    }
    return r;
}

std::vector<Pt> Shape::outline() const {
    if (kind == ShapeKind::polygon) return points;
    const Pt lo = points[0], hi = points[1];
    return {lo, {hi.x, lo.y}, hi, {lo.x, hi.y}};
}

int Geometry::background() const noexcept {
    for (std::size_t k = 0; k < regions.size(); ++k) {
        if (regions[k].shapes.empty()) return static_cast<int>(k);
    }
    return -1;
}

std::string format_length(Coord q) {
    if (q == 0) return "0";
    const bool negative = q < 0;
    const auto a = static_cast<std::uint64_t>(negative ? -q : q);
    std::string s = (negative ? "-" : "") + std::to_string(a / per_um);
    if (std::uint64_t frac = a % per_um; frac != 0) {
        std::string f = std::format("{:04}", frac);
        while (f.back() == '0') f.pop_back();
        s += "." + f;
    }
    return s + "um";
}

std::optional<Coord> parse_length(std::string_view s) {
    s = trim(s);
    double scale = 1e8;  // cm to quanta
    if (s.ends_with("um")) {
        scale = 1e4;
        s.remove_suffix(2);
    } else if (s.ends_with("nm")) {
        scale = 10.0;
        s.remove_suffix(2);
    } else if (s.ends_with("cm")) {
        s.remove_suffix(2);
    }
    if (!s.empty() && s.front() == '+') s.remove_prefix(1);
    double v = 0.0;
    const auto [end, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if (s.empty() || ec != std::errc{} || end != s.data() + s.size() || !std::isfinite(v)) return std::nullopt;
    return static_cast<Coord>(std::llround(v * scale));
}

// --- reading ---------------------------------------------------------------------------------------------------

std::expected<Geometry, std::string> read_geometry(std::string_view text, std::optional<Rect> bounds) {
    const auto d = deck::parse(text);
    if (!d) return std::unexpected(d.error().message);
    Geometry g;
    const deck::Section* mesh = nullptr;
    for (const auto& s : d->sections) {
        if (s.kind == "mesh") mesh = &s;
    }
    const auto refuse = [&](std::string why) {
        g.mesh = MeshKind::unsupported;
        g.why = std::move(why);
        return g;
    };
    if (!mesh) return refuse("the deck has no [mesh] section");
    const deck::Entry* type = find_entry(*mesh, "type");
    const std::string kind = type ? std::string(trim(type->value)) : std::string("tensor");
    if (kind == "tensor") {
        if (find_entry(*mesh, "z")) return refuse("3D geometry editing is a later unit");
        if (!find_entry(*mesh, "y")) return refuse("a 1D device has no geometry to draw");
        g.mesh = MeshKind::tensor;
    } else if (kind == "generate2d") {
        const deck::Entry* l = find_entry(*mesh, "layered");
        g.mesh = l && on_off(l->value).value_or(false) ? MeshKind::layered : MeshKind::tiled;
    } else {
        if (kind == "generate3d" || kind == "extrude3d") return refuse("3D geometry editing is a later unit");
        return refuse("a " + kind + " mesh is not drawn in the editor (2D tensor and generate2d decks are)");
    }
    std::map<std::string, std::pair<Pt, Pt>> patches;  // generate2d's named segments
    for (const auto& e : mesh->entries) {
        if (e.key.starts_with("patch ")) {
            const auto p = parse_points(e.value);
            if (!p || p->size() != 2) return refuse("line " + std::to_string(e.line) + ": a patch the editor cannot read");
            patches[e.key.substr(6)] = {(*p)[0], (*p)[1]};
        } else if (e.key == "domain") {
            const auto r = parse_box(e.value);
            if (!r) return refuse("line " + std::to_string(e.line) + ": a domain the editor cannot read");
            g.domain = *r;
            g.domain_given = true;
        }
    }
    int region_index = 0, contact_index = 0;
    for (const auto& s : d->sections) {
        if (s.kind == "region") {
            Region r{s.name, {}, {}, region_index++};
            for (const auto& e : s.entries) {
                if (e.key == "material") r.material = std::string(trim(e.value));
                if (!is_shape_key(e.key)) continue;
                auto shape = parse_shape(e);
                if (!shape) return refuse("line " + std::to_string(e.line) + ": a shape the editor cannot read");
                r.shapes.push_back(std::move(*shape));
            }
            g.regions.push_back(std::move(r));
        }
    }
    if (!g.domain_given) {
        if (bounds) {
            g.domain = *bounds;
        } else {
            bool any = false;
            for (const Region& r : g.regions) {
                for (const Shape& s : r.shapes) {
                    const Rect b = s.bounds();
                    g.domain = any ? normalized({std::min(g.domain.lo.x, b.lo.x), std::min(g.domain.lo.y, b.lo.y)},
                                                {std::max(g.domain.hi.x, b.hi.x), std::max(g.domain.hi.y, b.hi.y)})
                                   : b;
                    any = true;
                }
            }
        }
    }
    for (const auto& s : d->sections) {
        if (s.kind != "contact") continue;
        Contact c{s.name, {}, {}, {}, true, contact_index++};
        if (const auto* k = find_entry(s, "kind")) c.kind = std::string(trim(k->value));
        const deck::Entry* patch = find_entry(s, "patch");
        const std::string pname = patch ? std::string(trim(patch->value)) : std::string();
        if (g.mesh == MeshKind::tensor) {
            const auto face = face_segment(g.domain, pname);
            if (!face) {
                c.editable = false;
            } else {
                c.a = face->first;
                c.b = face->second;
                if (const auto* b = find_entry(s, "box")) {
                    const auto r = parse_box(b->value);
                    if (!r) {
                        c.editable = false;
                    } else if (c.a.x == c.b.x) {
                        c.a.y = std::max(c.a.y, r->lo.y);
                        c.b.y = std::min(c.b.y, r->hi.y);
                    } else {
                        c.a.x = std::max(c.a.x, r->lo.x);
                        c.b.x = std::min(c.b.x, r->hi.x);
                    }
                }
            }
        } else if (const auto it = patches.find(pname); it != patches.end()) {
            c.a = it->second.first;
            c.b = it->second.second;
        } else {
            c.editable = false;  // a face patch of another form
        }
        g.contacts.push_back(std::move(c));
    }
    return g;
}

// --- writing ---------------------------------------------------------------------------------------------------

std::expected<std::string, std::string> write_geometry(std::string_view text, const Geometry& g) {
    if (!g.editable()) return std::unexpected("this deck's geometry is not editable: " + g.why);
    const auto d = deck::parse(text);
    if (!d) return std::unexpected(d.error().message);
    const Lines l = split(text);
    std::vector<Change> changes;
    const auto block_of = [&](const deck::Section& s) {
        const int first = s.line - 1, last = section_last(l, s);
        return std::vector<std::string>(l.body.begin() + first, l.body.begin() + last + 1);
    };
    std::vector<const deck::Section*> regions, contacts;
    const deck::Section* mesh = nullptr;
    for (const auto& s : d->sections) {
        if (s.kind == "region") regions.push_back(&s);
        if (s.kind == "contact") contacts.push_back(&s);
        if (s.kind == "mesh") mesh = &s;
    }
    if (mesh == nullptr) return std::unexpected("the deck has no [mesh] section");

    // A section's lines with entries replaced (`replace`: the new text or nullopt to drop, per entry index), new
    // lines inserted after an entry (`after`, -1: after the header), and the header replaced.
    struct Edit {
        std::map<std::size_t, std::optional<std::string>> replace;
        std::multimap<int, std::string> after;
        std::optional<std::string> header;
    };
    const auto rewrite = [&](const deck::Section& s, const Edit& e) {
        std::vector<std::string> out;
        out.push_back(e.header ? *e.header : l.body[static_cast<std::size_t>(s.line - 1)]);
        const auto put_after = [&](int k) {
            for (auto [it, end] = e.after.equal_range(k); it != end; ++it) out.push_back(it->second);
        };
        put_after(-1);
        int line = s.line;  // 0-based index of the next line to copy
        for (std::size_t k = 0; k < s.entries.size(); ++k) {
            const int first = s.entries[k].line - 1, last = entry_last(l, first);
            for (; line < first; ++line) out.push_back(l.body[static_cast<std::size_t>(line)]);  // blanks, comments
            if (const auto it = e.replace.find(k); it != e.replace.end()) {
                if (it->second) out.push_back(*it->second);
            } else {
                for (int j = first; j <= last; ++j) out.push_back(l.body[static_cast<std::size_t>(j)]);
            }
            line = last + 1;
            put_after(static_cast<int>(k));
        }
        return out;
    };
    const auto change = [&](const deck::Section& s, std::vector<std::string> lines) {
        const int first = s.line - 1, last = section_last(l, s);
        if (lines != block_of(s)) changes.push_back({first, last - first + 1, std::move(lines)});
    };
    // A section removed with the blank line before it.
    const auto drop = [&](const deck::Section& s) {
        int first = s.line - 1;
        const int last = section_last(l, s);
        if (first > 0 && trim(l.body[static_cast<std::size_t>(first - 1)]).empty()) --first;
        changes.push_back({first, last - first + 1, {}});
    };

    // Regions: slot k (the k-th region section's place) takes the k-th region.
    const auto region_block = [&](const Region& r) -> std::vector<std::string> {
        if (r.origin < 0 || r.origin >= static_cast<int>(regions.size())) {
            std::vector<std::string> out{"", "[region " + r.name + "]"};
            if (!r.material.empty()) out.push_back("material = " + r.material);
            for (const Shape& s : r.shapes) out.push_back(shape_line(s));
            return out;
        }
        const deck::Section& s = *regions[static_cast<std::size_t>(r.origin)];
        Edit e;
        if (r.name != s.name) e.header = "[region " + r.name + "]";
        std::size_t shape = 0;
        int last_shape = -1, last_entry = s.entries.empty() ? -1 : static_cast<int>(s.entries.size()) - 1;
        bool material = false;
        for (std::size_t k = 0; k < s.entries.size(); ++k) {
            const auto& en = s.entries[k];
            if (en.key == "material") {
                material = true;
                if (std::string(trim(en.value)) != r.material) {
                    e.replace[k] = r.material.empty() ? std::nullopt : std::optional<std::string>("material = " + r.material);
                }
            } else if (is_shape_key(en.key)) {
                last_shape = static_cast<int>(k);
                if (shape >= r.shapes.size()) {
                    e.replace[k] = std::nullopt;
                } else if (parse_shape(en) != r.shapes[shape]) {
                    e.replace[k] = shape_line(r.shapes[shape]);
                }
                ++shape;
            }
        }
        if (!material && !r.material.empty()) e.after.emplace(-1, "material = " + r.material);
        for (; shape < r.shapes.size(); ++shape) e.after.emplace(last_shape >= 0 ? last_shape : last_entry, shape_line(r.shapes[shape]));
        return rewrite(s, e);
    };
    for (std::size_t k = 0; k < std::max(regions.size(), g.regions.size()); ++k) {
        if (k < regions.size() && k < g.regions.size()) {
            auto block = region_block(g.regions[k]);
            if (g.regions[k].origin < 0) block.erase(block.begin());  // a slot keeps its own spacing
            change(*regions[k], std::move(block));
        } else if (k < regions.size()) {
            drop(*regions[k]);
        } else {
            const int at = regions.empty() ? section_last(l, *mesh) + 1 : section_last(l, *regions.back()) + 1;
            changes.push_back({at, 0, region_block(Region{g.regions[k].name, g.regions[k].material, g.regions[k].shapes, -1})});
        }
    }

    // Each contact's patch name: a kept contact's as its section says; a new one's its name, made unique among the
    // deck's patches.
    std::vector<std::string> patch_of(g.contacts.size());
    {
        std::set<std::string> taken;
        for (const auto& en : mesh->entries) {
            if (en.key.starts_with("patch ")) taken.insert(en.key.substr(6));
        }
        for (std::size_t i = 0; i < g.contacts.size(); ++i) {
            const Contact& c = g.contacts[i];
            if (c.origin >= 0 && c.origin < static_cast<int>(contacts.size())) {
                if (const auto* p = find_entry(*contacts[static_cast<std::size_t>(c.origin)], "patch")) patch_of[i] = std::string(trim(p->value));
                taken.insert(patch_of[i]);
            }
        }
        for (std::size_t i = 0; i < g.contacts.size(); ++i) {
            if (!patch_of[i].empty()) continue;
            std::string name = g.contacts[i].name;
            for (int n = 2; taken.contains(name); ++n) name = g.contacts[i].name + "_" + std::to_string(n);
            patch_of[i] = name;
            taken.insert(name);
        }
    }

    // [mesh]: layered, domain, the generated deck's patches.
    {
        Edit e;
        int type_entry = -1, domain_entry = -1, last_patch = -1;
        bool layered = false;
        std::map<std::string, std::size_t> patch_entries;
        for (std::size_t k = 0; k < mesh->entries.size(); ++k) {
            const auto& en = mesh->entries[k];
            if (en.key == "type") type_entry = static_cast<int>(k);
            if (en.key == "layered") {
                layered = true;
                if (g.mesh == MeshKind::layered && !on_off(en.value).value_or(false)) e.replace[k] = "layered = on";
            }
            if (en.key == "domain") domain_entry = static_cast<int>(k);
            if (en.key.starts_with("patch ")) {
                patch_entries[en.key.substr(6)] = k;
                last_patch = static_cast<int>(k);
            }
        }
        if (g.mesh == MeshKind::layered && !layered) e.after.emplace(type_entry, "layered = on");
        const std::string domain_line = "domain = " + format_length(g.domain.lo.x) + " " + format_length(g.domain.hi.x) + " " +
                                        format_length(g.domain.lo.y) + " " + format_length(g.domain.hi.y);
        if (domain_entry >= 0) {
            const auto& en = mesh->entries[static_cast<std::size_t>(domain_entry)];
            if (!g.domain_given) e.replace[static_cast<std::size_t>(domain_entry)] = std::nullopt;
            else if (parse_box(en.value) != g.domain) e.replace[static_cast<std::size_t>(domain_entry)] = domain_line;
        } else if (g.domain_given) {
            e.after.emplace(type_entry, domain_line);
        }
        if (g.mesh != MeshKind::tensor) {
            // Each kept contact's patch; a patch no kept contact uses is dropped; a new contact's patch added.
            std::set<std::string> used;
            for (std::size_t i = 0; i < g.contacts.size(); ++i) {
                const Contact& c = g.contacts[i];
                const std::string& pname = patch_of[i];
                if (!c.editable) {
                    used.insert(pname);
                    continue;
                }
                const std::string line = "patch " + pname + " = " + segment_text(c.a, c.b);
                if (const auto it = patch_entries.find(pname); it != patch_entries.end()) {
                    const auto p = parse_points(mesh->entries[it->second].value);
                    if (!used.contains(pname) && (!p || p->size() != 2 || (*p)[0] != c.a || (*p)[1] != c.b)) e.replace[it->second] = line;
                } else {
                    e.after.emplace(last_patch >= 0 ? last_patch : type_entry, line);
                }
                used.insert(pname);
            }
            for (const auto& [name, k] : patch_entries) {
                if (!used.contains(name)) e.replace[k] = std::nullopt;
            }
        }
        change(*mesh, rewrite(*mesh, e));
    }

    // Contacts: kept ones renamed or moved in place; removed ones dropped; new ones appended.
    std::vector<bool> kept(contacts.size(), false);
    std::vector<std::string> appended;
    for (const Contact& c : g.contacts) {
        if (c.origin < 0 || c.origin >= static_cast<int>(contacts.size())) {
            appended.insert(appended.end(), {"", "[contact " + c.name + "]", "kind = " + (c.kind.empty() ? std::string("ohmic") : c.kind)});
            if (g.mesh == MeshKind::tensor) {
                appended.push_back("patch = " + face_of(g.domain, c.a, c.b));
                appended.push_back("box = " + format_length(std::min(c.a.x, c.b.x)) + " " + format_length(std::max(c.a.x, c.b.x)) + " " +
                                   format_length(std::min(c.a.y, c.b.y)) + " " + format_length(std::max(c.a.y, c.b.y)));
            } else {
                appended.push_back("patch = " + patch_of[static_cast<std::size_t>(&c - g.contacts.data())]);
            }
            continue;
        }
        kept[static_cast<std::size_t>(c.origin)] = true;
        const deck::Section& s = *contacts[static_cast<std::size_t>(c.origin)];
        Edit e;
        if (c.name != s.name) e.header = "[contact " + c.name + "]";
        for (std::size_t k = 0; k < s.entries.size(); ++k) {
            if (s.entries[k].key == "kind" && std::string(trim(s.entries[k].value)) != c.kind && !c.kind.empty()) {
                e.replace[k] = "kind = " + c.kind;
            }
        }
        if (g.mesh == MeshKind::tensor && c.editable) {
            const std::string face = face_of(g.domain, c.a, c.b);
            const auto whole = face_segment(g.domain, face);
            const Rect want = normalized(c.a, c.b);
            const std::string box_line = "box = " + format_length(want.lo.x) + " " + format_length(want.hi.x) + " " +
                                         format_length(want.lo.y) + " " + format_length(want.hi.y);
            int box_entry = -1, patch_entry = -1;
            for (std::size_t k = 0; k < s.entries.size(); ++k) {
                if (s.entries[k].key == "box") box_entry = static_cast<int>(k);
                if (s.entries[k].key == "patch") patch_entry = static_cast<int>(k);
            }
            if (patch_entry >= 0 && std::string(trim(s.entries[static_cast<std::size_t>(patch_entry)].value)) != face) {
                e.replace[static_cast<std::size_t>(patch_entry)] = "patch = " + face;
            }
            // The box as read back: the face clipped to it must give the segment.
            const bool whole_face = whole && normalized(whole->first, whole->second) == want;
            if (box_entry >= 0) {
                const auto r = parse_box(s.entries[static_cast<std::size_t>(box_entry)].value);
                Contact again = c;
                if (r && whole) {
                    again.a = whole->first;
                    again.b = whole->second;
                    if (again.a.x == again.b.x) {
                        again.a.y = std::max(again.a.y, r->lo.y);
                        again.b.y = std::min(again.b.y, r->hi.y);
                    } else {
                        again.a.x = std::max(again.a.x, r->lo.x);
                        again.b.x = std::min(again.b.x, r->hi.x);
                    }
                }
                if (!r || normalized(again.a, again.b) != want) e.replace[static_cast<std::size_t>(box_entry)] = box_line;
            } else if (!whole_face) {
                e.after.emplace(patch_entry, box_line);
            }
        }
        change(s, rewrite(s, e));
    }
    for (std::size_t k = 0; k < contacts.size(); ++k) {
        if (!kept[k]) {
            drop(*contacts[k]);
        }
    }
    if (!appended.empty()) changes.push_back({static_cast<int>(l.body.size()), 0, std::move(appended)});

    // Applied in order; an insertion keeps its order among those at one place.
    std::stable_sort(changes.begin(), changes.end(), [](const Change& a, const Change& b) { return a.first < b.first; });
    std::string out;
    out.reserve(text.size() + 256);
    int line = 0;
    const int n = static_cast<int>(l.body.size());
    const auto copy = [&](int k) { out += l.body[static_cast<std::size_t>(k)] + l.eol[static_cast<std::size_t>(k)]; };
    for (const Change& c : changes) {
        for (; line < c.first && line < n; ++line) copy(line);
        if (c.first >= n && !out.empty() && !out.ends_with('\n')) out += l.nl;  // after a last line without one
        for (const std::string& s : c.lines) out += s + l.nl;
        line = std::max(line, c.first + c.count);
    }
    for (; line < n; ++line) copy(line);
    return out;
}

// --- edits -----------------------------------------------------------------------------------------------------

namespace {

// A polygon the layered kernel takes (slice 4, GD10: edges in any direction): three or more vertices within ±10 cm,
// no edge empty, simple, of non-zero area. Exact: with coordinates within 1e9 quanta every cross product below fits
// in 64 bits.
constexpr Coord coord_limit = 1'000'000'000;
int turn(Pt a, Pt b, Pt c) {
    const std::int64_t v = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    return v > 0 ? 1 : v < 0 ? -1 : 0;
}
bool on_closed(Pt a, Pt b, Pt p) {
    return turn(a, b, p) == 0 && p.x >= std::min(a.x, b.x) && p.x <= std::max(a.x, b.x) && p.y >= std::min(a.y, b.y) &&
           p.y <= std::max(a.y, b.y);
}
bool meet(Pt a, Pt b, Pt c, Pt d) {
    const int o1 = turn(a, b, c), o2 = turn(a, b, d), o3 = turn(c, d, a), o4 = turn(c, d, b);
    if (o1 * o2 < 0 && o3 * o4 < 0) return true;
    return on_closed(a, b, c) || on_closed(a, b, d) || on_closed(c, d, a) || on_closed(c, d, b);
}
bool outline_ok(const std::vector<Pt>& p) {
    const std::size_t n = p.size();
    if (n < 3) return false;
    for (std::size_t i = 0; i < n; ++i) {
        const Pt a = p[i], b = p[(i + 1) % n], c = p[(i + 2) % n];
        if (std::abs(a.x) > coord_limit || std::abs(a.y) > coord_limit || a == b) return false;
        // Consecutive edges meet only at their shared vertex (no reversal: an outline without area turns back
        // somewhere, so this refuses it too); the others not at all.
        if (turn(a, b, c) == 0 && on_closed(a, b, c)) return false;
        for (std::size_t j = i + 2; j < n; ++j) {
            if (i == 0 && j == n - 1) continue;
            if (meet(a, b, p[j], p[(j + 1) % n])) return false;
        }
    }
    return true;
}

// Whether `g` can be edited; a tiled generated deck becomes layered (its polygons tile, so painting them in turn
// gives the same device; slanted edges too since slice 4).
bool editing(Geometry& g) {
    if (g.mesh == MeshKind::tiled) g.mesh = MeshKind::layered;
    return g.editable();
}

bool valid(Geometry& g, ShapeRef s) {
    return editing(g) && s.region >= 0 && s.region < static_cast<int>(g.regions.size()) && s.shape >= 0 &&
           s.shape < static_cast<int>(g.regions[static_cast<std::size_t>(s.region)].shapes.size());
}

Shape& at(Geometry& g, ShapeRef s) {
    return g.regions[static_cast<std::size_t>(s.region)].shapes[static_cast<std::size_t>(s.shape)];
}

bool identifier(std::string_view s) {
    if (s.empty() || !std::isalpha(static_cast<unsigned char>(s.front()))) return false;
    return std::ranges::all_of(s, [](char c) { return std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '.' || c == '-'; });
}

std::string fresh(const Geometry& g, std::string base) {
    if (!identifier(base)) base = "region";
    for (int n = 1;; ++n) {
        std::string name = base + "_" + std::to_string(n);
        if (std::ranges::none_of(g.regions, [&](const Region& r) { return r.name == name; })) return name;
    }
}

int add_region(Geometry& g, Shape s, std::string material) {
    Region r{fresh(g, material), std::move(material), {std::move(s)}, -1};
    g.regions.push_back(std::move(r));
    return static_cast<int>(g.regions.size()) - 1;
}

}  // namespace

int add_box(Geometry& g, Rect r, std::string material) {
    if (!editing(g) || r.lo.x >= r.hi.x || r.lo.y >= r.hi.y) return -1;
    return add_region(g, Shape{ShapeKind::box, {r.lo, r.hi}}, std::move(material));
}

int add_polygon(Geometry& g, std::vector<Pt> points, std::string material) {
    if (!editing(g) || g.mesh == MeshKind::tensor || !outline_ok(points)) return -1;
    return add_region(g, Shape{ShapeKind::polygon, std::move(points)}, std::move(material));
}

bool move_shape(Geometry& g, ShapeRef s, Coord dx, Coord dy) {
    if (!valid(g, s)) return false;
    for (Pt& p : at(g, s).points) p = {p.x + dx, p.y + dy};
    return true;
}

bool set_box(Geometry& g, ShapeRef s, Rect r) {
    if (!valid(g, s) || at(g, s).kind != ShapeKind::box || r.lo.x >= r.hi.x || r.lo.y >= r.hi.y) return false;
    at(g, s).points = {r.lo, r.hi};
    return true;
}

bool set_vertex(Geometry& g, ShapeRef s, int k, Pt p, bool follow) {
    if (!valid(g, s) || at(g, s).kind != ShapeKind::polygon) return false;
    auto& v = at(g, s).points;
    if (k < 0 || k >= static_cast<int>(v.size())) return false;
    // With `follow`, its neighbours follow along the horizontal and vertical edges they share with it, so those stay
    // so; a slanted edge's other end stays.
    auto w = v;
    const std::size_t n = w.size(), i = static_cast<std::size_t>(k);
    const Pt old = w[i];
    if (follow) {
        for (const std::size_t j : {(i + n - 1) % n, (i + 1) % n}) {
            if (w[j].x == old.x) w[j].x = p.x;
            else if (w[j].y == old.y) w[j].y = p.y;
        }
    }
    w[i] = p;
    if (!outline_ok(w)) return false;
    v = std::move(w);
    return true;
}

bool insert_vertex(Geometry& g, ShapeRef s, int after, Pt p) {  // on the edge after `after`
    if (!valid(g, s) || at(g, s).kind != ShapeKind::polygon) return false;
    auto& v = at(g, s).points;
    if (after < 0 || after >= static_cast<int>(v.size())) return false;
    const Pt a = v[static_cast<std::size_t>(after)], b = v[(static_cast<std::size_t>(after) + 1) % v.size()];
    if (!on_closed(a, b, p) || p == a || p == b) return false;
    v.insert(v.begin() + after + 1, p);
    return true;
}

bool remove_vertex(Geometry& g, ShapeRef s, int k) {  // the polygon left valid (three vertices or more)
    if (!valid(g, s) || at(g, s).kind != ShapeKind::polygon) return false;
    auto& v = at(g, s).points;
    if (k < 0 || k >= static_cast<int>(v.size())) return false;
    auto w = v;
    w.erase(w.begin() + k);
    if (!outline_ok(w)) return false;
    v = std::move(w);
    return true;
}

bool delete_shape(Geometry& g, ShapeRef s) {
    if (!valid(g, s)) return false;
    auto& shapes = g.regions[static_cast<std::size_t>(s.region)].shapes;
    shapes.erase(shapes.begin() + s.shape);
    if (shapes.empty()) g.regions.erase(g.regions.begin() + s.region);
    return true;
}

bool rename_region(Geometry& g, int region, std::string name) {
    if (!editing(g) || region < 0 || region >= static_cast<int>(g.regions.size()) || !identifier(name)) return false;
    for (std::size_t k = 0; k < g.regions.size(); ++k) {
        if (static_cast<int>(k) != region && g.regions[k].name == name) return false;
    }
    g.regions[static_cast<std::size_t>(region)].name = std::move(name);
    return true;
}

bool set_material(Geometry& g, int region, std::string material) {
    if (!editing(g) || region < 0 || region >= static_cast<int>(g.regions.size()) || !identifier(material)) return false;
    g.regions[static_cast<std::size_t>(region)].material = std::move(material);
    return true;
}

bool raise_region(Geometry& g, int region) {
    if (!editing(g) || region < 0 || region + 1 >= static_cast<int>(g.regions.size())) return false;
    std::swap(g.regions[static_cast<std::size_t>(region)], g.regions[static_cast<std::size_t>(region + 1)]);
    return true;
}

bool lower_region(Geometry& g, int region) {
    if (!editing(g) || region <= 0 || region >= static_cast<int>(g.regions.size())) return false;
    return raise_region(g, region - 1);
}

int duplicate_shape(Geometry& g, ShapeRef s, Coord dx, Coord dy) {
    if (!valid(g, s)) return -1;
    Shape copy = at(g, s);
    for (Pt& p : copy.points) p = {p.x + dx, p.y + dy};
    return add_region(g, std::move(copy), g.regions[static_cast<std::size_t>(s.region)].material);
}

bool mirror_shape(Geometry& g, ShapeRef s, bool horizontal) {
    if (!valid(g, s)) return false;
    Shape& sh = at(g, s);
    const Rect b = sh.bounds();
    for (Pt& p : sh.points) {
        if (horizontal) p.x = b.lo.x + b.hi.x - p.x;
        else p.y = b.lo.y + b.hi.y - p.y;
    }
    if (sh.kind == ShapeKind::box) {
        const Rect r = normalized(sh.points[0], sh.points[1]);
        sh.points = {r.lo, r.hi};
    } else {
        std::ranges::reverse(sh.points);  // keep the orientation
    }
    return true;
}

bool on_one_side(const Rect& d, Pt a, Pt b) noexcept {
    if (a == b) return false;
    if (a.x == b.x && (a.x == d.lo.x || a.x == d.hi.x)) return std::min(a.y, b.y) >= d.lo.y && std::max(a.y, b.y) <= d.hi.y;
    if (a.y == b.y && (a.y == d.lo.y || a.y == d.hi.y)) return std::min(a.x, b.x) >= d.lo.x && std::max(a.x, b.x) <= d.hi.x;
    return false;
}

int add_contact(Geometry& g, std::string name, std::string kind, Pt a, Pt b) {
    if (!editing(g) || !identifier(name) || !on_one_side(g.domain, a, b)) return -1;
    if (std::ranges::any_of(g.contacts, [&](const Contact& c) { return c.name == name; })) return -1;
    if (b < a) std::swap(a, b);
    g.contacts.push_back({std::move(name), std::move(kind), a, b, true, -1});
    return static_cast<int>(g.contacts.size()) - 1;
}

bool set_contact(Geometry& g, int contact, Pt a, Pt b) {
    if (!editing(g) || contact < 0 || contact >= static_cast<int>(g.contacts.size())) return false;
    Contact& c = g.contacts[static_cast<std::size_t>(contact)];
    if (!c.editable || !on_one_side(g.domain, a, b)) return false;
    if (b < a) std::swap(a, b);
    c.a = a;
    c.b = b;
    return true;
}

bool delete_contact(Geometry& g, int contact) {
    if (!editing(g) || contact < 0 || contact >= static_cast<int>(g.contacts.size())) return false;
    g.contacts.erase(g.contacts.begin() + contact);
    return true;
}

bool rename_contact(Geometry& g, int contact, std::string name) {
    if (!editing(g) || contact < 0 || contact >= static_cast<int>(g.contacts.size()) || !identifier(name)) return false;
    for (std::size_t k = 0; k < g.contacts.size(); ++k) {
        if (static_cast<int>(k) != contact && g.contacts[k].name == name) return false;
    }
    g.contacts[static_cast<std::size_t>(contact)].name = std::move(name);
    return true;
}

// --- undo ------------------------------------------------------------------------------------------------------

void EditHistory::record(std::string before) {
    undo_.push_back(std::move(before));
    if (undo_.size() > limit) undo_.erase(undo_.begin());
    redo_.clear();
}

std::optional<std::string> EditHistory::undo(std::string current) {
    if (undo_.empty()) return std::nullopt;
    std::string t = std::move(undo_.back());
    undo_.pop_back();
    redo_.push_back(std::move(current));
    return t;
}

std::optional<std::string> EditHistory::redo(std::string current) {
    if (redo_.empty()) return std::nullopt;
    std::string t = std::move(redo_.back());
    redo_.pop_back();
    undo_.push_back(std::move(current));
    return t;
}

void EditHistory::clear() noexcept {
    undo_.clear();
    redo_.clear();
}

// --- view, snapping, hit tests ---------------------------------------------------------------------------------

Pt View::model(double x, double y) const noexcept {
    return {static_cast<Coord>(std::llround((x - origin_x) / px_per_um * per_um)),
            static_cast<Coord>(std::llround((origin_y - y) / px_per_um * per_um))};
}

Coord grid_step(const View& v) noexcept {
    const double min_px = v.px(12.0);
    for (Coord decade = 1;; decade *= 10) {
        for (const Coord m : {1, 2, 5}) {
            const Coord step = m * decade;
            if (static_cast<double>(step) / per_um * v.px_per_um >= min_px) return step;
        }
        if (decade > Coord{1} << 50) return decade;
    }
}

namespace {

struct Seg {
    Pt a, b;
};

// The edges of every shape (but `skip`) and the domain's sides.
std::vector<Seg> edges(const Geometry& g, std::optional<ShapeRef> skip) {
    std::vector<Seg> out;
    for (std::size_t r = 0; r < g.regions.size(); ++r) {
        for (std::size_t s = 0; s < g.regions[r].shapes.size(); ++s) {
            if (skip && skip->region == static_cast<int>(r) && skip->shape == static_cast<int>(s)) continue;
            const auto o = g.regions[r].shapes[s].outline();
            for (std::size_t k = 0; k < o.size(); ++k) out.push_back({o[k], o[(k + 1) % o.size()]});
        }
    }
    const auto d = Shape{ShapeKind::box, {g.domain.lo, g.domain.hi}}.outline();
    for (std::size_t k = 0; k < 4; ++k) out.push_back({d[k], d[(k + 1) % 4]});
    return out;
}

double dist2(double ax, double ay, double bx, double by) { return (ax - bx) * (ax - bx) + (ay - by) * (ay - by); }

// The nearest point of segment ab to p, in screen space; t its parameter.
std::pair<double, double> nearest(const View& v, Seg s, double x, double y, double* t = nullptr) {
    const double ax = v.sx(s.a.x), ay = v.sy(s.a.y), bx = v.sx(s.b.x), by = v.sy(s.b.y);
    const double dx = bx - ax, dy = by - ay, len2 = dx * dx + dy * dy;
    double u = len2 > 0.0 ? ((x - ax) * dx + (y - ay) * dy) / len2 : 0.0;
    u = std::clamp(u, 0.0, 1.0);
    if (t) *t = u;
    return {ax + u * dx, ay + u * dy};
}

bool inside(const std::vector<Pt>& poly, double x, double y) {  // model doubles, even-odd
    bool in = false;
    for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
        const double xi = static_cast<double>(poly[i].x), yi = static_cast<double>(poly[i].y);
        const double xj = static_cast<double>(poly[j].x), yj = static_cast<double>(poly[j].y);
        if ((yi > y) != (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) in = !in;
    }
    return in;
}

}  // namespace

Snap snap(const Geometry& g, const View& v, Pt raw, std::optional<ShapeRef> skip, bool grid) {
    const double tol2 = v.px(8.0) * v.px(8.0);
    const double x = v.sx(raw.x), y = v.sy(raw.y);
    const auto segs = edges(g, skip);
    // Vertices: the shapes', the contacts' ends, the domain's corners.
    std::optional<Snap> best;
    double best_d = tol2;
    const auto consider = [&](Pt p, SnapKind k) {
        const double d = dist2(v.sx(p.x), v.sy(p.y), x, y);
        if (d <= best_d) {
            best_d = d;
            best = Snap{p, k};
        }
    };
    for (const Seg& s : segs) consider(s.a, SnapKind::vertex);
    for (const Contact& c : g.contacts) {
        consider(c.a, SnapKind::vertex);
        consider(c.b, SnapKind::vertex);
    }
    if (best) return *best;
    for (const Seg& s : segs) consider({(s.a.x + s.b.x) / 2, (s.a.y + s.b.y) / 2}, SnapKind::midpoint);
    if (best) return *best;
    for (const Seg& s : segs) {
        const auto [px, py] = nearest(v, s, x, y);
        const double d = dist2(px, py, x, y);
        if (d <= best_d) {
            best_d = d;
            Pt p = v.model(px, py);
            if (s.a.x == s.b.x) p.x = s.a.x;  // exactly on an axis-parallel edge
            if (s.a.y == s.b.y) p.y = s.a.y;
            best = Snap{p, SnapKind::edge};
        }
    }
    if (best) return *best;
    if (!grid) return {raw, SnapKind::none};
    const Coord step = grid_step(v);
    const auto round = [&](Coord c) { return static_cast<Coord>(std::llround(static_cast<double>(c) / static_cast<double>(step))) * step; };
    return {{round(raw.x), round(raw.y)}, SnapKind::grid};
}

Pt constrain(Pt anchor, Pt p) noexcept {
    const Coord dx = p.x - anchor.x, dy = p.y - anchor.y;
    return (dx < 0 ? -dx : dx) >= (dy < 0 ? -dy : dy) ? Pt{p.x, anchor.y} : Pt{anchor.x, p.y};
}

Hit hit_test(const Geometry& g, const View& v, double x, double y, std::optional<ShapeRef> selected) {
    if (selected && selected->region >= 0 && selected->region < static_cast<int>(g.regions.size()) && selected->shape >= 0 &&
        selected->shape < static_cast<int>(g.regions[static_cast<std::size_t>(selected->region)].shapes.size())) {
        const Shape& s = g.regions[static_cast<std::size_t>(selected->region)].shapes[static_cast<std::size_t>(selected->shape)];
        const double h = v.px(4.5);
        if (s.kind == ShapeKind::box) {
            const Rect r = normalized(s.points[0], s.points[1]);
            const Coord mx = (r.lo.x + r.hi.x) / 2, my = (r.lo.y + r.hi.y) / 2;
            const Pt handles[8] = {r.lo, {mx, r.lo.y}, {r.hi.x, r.lo.y}, {r.hi.x, my}, r.hi, {mx, r.hi.y}, {r.lo.x, r.hi.y}, {r.lo.x, my}};
            for (int k = 0; k < 8; ++k) {
                if (std::abs(v.sx(handles[k].x) - x) <= h && std::abs(v.sy(handles[k].y) - y) <= h) return {HitKind::handle, *selected, k, -1};
            }
        } else {
            for (std::size_t k = 0; k < s.points.size(); ++k) {
                if (std::abs(v.sx(s.points[k].x) - x) <= h && std::abs(v.sy(s.points[k].y) - y) <= h) {
                    return {HitKind::vertex, *selected, static_cast<int>(k), -1};
                }
            }
        }
    }
    const double c2 = v.px(5.0) * v.px(5.0);
    for (std::size_t k = 0; k < g.contacts.size(); ++k) {
        const Contact& c = g.contacts[k];
        for (int e = 0; e < 2; ++e) {
            const Pt p = e == 0 ? c.a : c.b;
            if (dist2(v.sx(p.x), v.sy(p.y), x, y) <= c2) return {HitKind::contact_end, {}, static_cast<int>(k), e};
        }
    }
    for (std::size_t k = 0; k < g.contacts.size(); ++k) {
        const auto [px, py] = nearest(v, {g.contacts[k].a, g.contacts[k].b}, x, y);
        if (dist2(px, py, x, y) <= c2) return {HitKind::contact, {}, static_cast<int>(k), -1};
    }
    if (selected && selected->region >= 0 && selected->region < static_cast<int>(g.regions.size()) && selected->shape >= 0 &&
        selected->shape < static_cast<int>(g.regions[static_cast<std::size_t>(selected->region)].shapes.size())) {
        const Shape& s = g.regions[static_cast<std::size_t>(selected->region)].shapes[static_cast<std::size_t>(selected->shape)];
        if (s.kind == ShapeKind::polygon) {
            const double e2 = v.px(4.0) * v.px(4.0);
            for (std::size_t k = 0; k < s.points.size(); ++k) {
                const auto [px, py] = nearest(v, {s.points[k], s.points[(k + 1) % s.points.size()]}, x, y);
                if (dist2(px, py, x, y) <= e2) return {HitKind::edge, *selected, static_cast<int>(k), -1};
            }
        }
    }
    // The topmost shape whose outline holds the point (model coordinates, unrounded).
    const double mx = (x - v.origin_x) / v.px_per_um * per_um, my = (v.origin_y - y) / v.px_per_um * per_um;
    for (std::size_t r = g.regions.size(); r-- > 0;) {
        const auto& shapes = g.regions[r].shapes;
        for (std::size_t s = shapes.size(); s-- > 0;) {
            if (inside(shapes[s].outline(), mx, my)) return {HitKind::shape, {static_cast<int>(r), static_cast<int>(s)}, -1, -1};
        }
    }
    const Rect& d = g.domain;
    if (const int bg = g.background(); bg >= 0 && mx >= static_cast<double>(d.lo.x) && mx <= static_cast<double>(d.hi.x) &&
                                       my >= static_cast<double>(d.lo.y) && my <= static_cast<double>(d.hi.y)) {
        return {HitKind::shape, {bg, -1}, -1, -1};
    }
    return {};
}

Rect drag_handle(const Rect& r, int h, Pt p) noexcept {
    Pt lo = r.lo, hi = r.hi;
    switch (h) {
    case 0: lo = p; break;
    case 1: lo.y = p.y; break;
    case 2: lo.y = p.y; hi.x = p.x; break;
    case 3: hi.x = p.x; break;
    case 4: hi = p; break;
    case 5: hi.y = p.y; break;
    case 6: lo.x = p.x; hi.y = p.y; break;
    case 7: lo.x = p.x; break;
    default: break;
    }
    return normalized(lo, hi);
}

}  // namespace NiTCAD::app::geo

namespace NiTCAD::app::geo {

bool set_contact_kind(Geometry& g, int contact, std::string kind) {
    if (!editing(g) || contact < 0 || contact >= static_cast<int>(g.contacts.size()) || !identifier(kind)) return false;
    g.contacts[static_cast<std::size_t>(contact)].kind = std::move(kind);
    return true;
}

bool set_domain(Geometry& g, Rect d) {
    if (!editing(g) || g.mesh != MeshKind::layered || d.lo.x >= d.hi.x || d.lo.y >= d.hi.y) return false;
    // Each contact moves with its side, clipped to it; refused when one would be empty.
    std::vector<Contact> moved = g.contacts;
    for (Contact& c : moved) {
        const Rect& o = g.domain;
        const bool vertical = c.a.x == c.b.x;
        if (vertical) {
            const Coord x = c.a.x == o.lo.x ? d.lo.x : d.hi.x;
            c.a = {x, std::max(c.a.y, d.lo.y)};
            c.b = {x, std::min(c.b.y, d.hi.y)};
        } else {
            const Coord y = c.a.y == o.lo.y ? d.lo.y : d.hi.y;
            c.a = {std::max(c.a.x, d.lo.x), y};
            c.b = {std::min(c.b.x, d.hi.x), y};
        }
        if (c.editable && !on_one_side(d, c.a, c.b)) return false;
    }
    g.contacts = std::move(moved);
    g.domain = d;
    g.domain_given = true;
    return true;
}

}  // namespace NiTCAD::app::geo

// --- validation (GD9) --------------------------------------------------------------------------------------------

namespace NiTCAD::app::geo {

namespace {

// The points "(x, y)" or "(x, y, z)" a message names, in cm (the mesher's and the build's form).
std::vector<Pt> named_points(std::string_view m) {
    std::vector<Pt> out;
    for (std::size_t open = m.find('('); open != std::string_view::npos; open = m.find('(', open + 1)) {
        const std::size_t close = m.find(')', open);
        if (close == std::string_view::npos) break;
        std::vector<double> v;
        std::string_view inner = m.substr(open + 1, close - open - 1);
        bool ok = true;
        while (ok && !inner.empty()) {
            const std::size_t comma = inner.find(',');
            const std::string_view t = trim(inner.substr(0, comma));
            double x = 0.0;
            const auto r = std::from_chars(t.data(), t.data() + t.size(), x);
            ok = !t.empty() && r.ec == std::errc{} && r.ptr == t.data() + t.size() && std::isfinite(x);
            if (ok) v.push_back(x);
            inner = comma == std::string_view::npos ? std::string_view{} : inner.substr(comma + 1);
        }
        if (ok && (v.size() == 2 || v.size() == 3)) {
            out.push_back({static_cast<Coord>(std::llround(v[0] / quantum_cm)), static_cast<Coord>(std::llround(v[1] / quantum_cm))});
        }
    }
    return out;
}

// Two contacts along one line of the boundary sharing more than an end point.
bool overlapping(const Contact& p, const Contact& q) noexcept {
    const bool pv = p.a.x == p.b.x, qv = q.a.x == q.b.x;
    if (pv != qv) return false;
    if (pv ? p.a.x != q.a.x : p.a.y != q.a.y) return false;
    const auto span = [&](const Contact& c) {
        return pv ? std::pair{std::min(c.a.y, c.b.y), std::max(c.a.y, c.b.y)} : std::pair{std::min(c.a.x, c.b.x), std::max(c.a.x, c.b.x)};
    };
    const auto [p0, p1] = span(p);
    const auto [q0, q1] = span(q);
    return std::max(p0, q0) < std::min(p1, q1) || (p0 == q0 && p1 == q1);
}

}  // namespace

Problem locate_problem(const Geometry& g, std::string_view text, int line, std::string message) {
    Problem out{line, std::move(message), {}};
    const auto d = deck::parse(text);
    if (d && line > 0) {
        // The section holding the line: the last one starting at or before it.
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
        const auto contact_of = [&](int origin) {
            for (std::size_t k = 0; k < g.contacts.size(); ++k) {
                if (g.contacts[k].origin == origin) return static_cast<int>(k);
            }
            return -1;
        };
        const auto add_contact = [&](int k) {
            if (k < 0) return;
            const Marker mk{Marker::Kind::contact, {}, k, {}};
            if (std::ranges::find(out.markers, mk) == out.markers.end()) out.markers.push_back(mk);
        };
        if (at && region_index >= 0) {
            for (std::size_t r = 0; r < g.regions.size(); ++r) {
                if (g.regions[r].origin != region_index) continue;
                int shape = -1, n = 0;
                for (const auto& e : at->entries) {
                    if (!is_shape_key(e.key)) continue;
                    if (e.line == line) shape = n;
                    ++n;
                }
                if (g.regions[r].shapes.empty()) out.markers.push_back({Marker::Kind::domain, {static_cast<int>(r), -1}, -1, {}});
                else out.markers.push_back({Marker::Kind::shape, {static_cast<int>(r), shape}, -1, {}});
            }
        } else if (at && contact_index >= 0) {
            add_contact(contact_of(contact_index));
        } else if (at && at->kind == "mesh") {
            const deck::Entry* patch = nullptr;
            for (const auto& e : at->entries) {
                if (e.line == line && e.key.starts_with("patch ")) patch = &e;
            }
            if (patch) {
                // The contacts using this patch, and every contact overlapping one of them.
                const std::string_view name = trim(std::string_view(patch->key).substr(6));
                int c = 0;
                for (const auto& s : d->sections) {
                    if (s.kind != "contact") continue;
                    const deck::Entry* p = find_entry(s, "patch");
                    if (p && trim(p->value) == name) add_contact(contact_of(c));
                    ++c;
                }
                const std::size_t direct = out.markers.size();
                for (std::size_t m = 0; m < direct; ++m) {
                    const Contact& p = g.contacts[static_cast<std::size_t>(out.markers[m].contact)];
                    for (std::size_t k = 0; k < g.contacts.size(); ++k) {
                        if (static_cast<int>(k) != out.markers[m].contact && overlapping(p, g.contacts[k])) add_contact(static_cast<int>(k));
                    }
                }
            } else if (out.message.find("(") == std::string::npos) {
                out.markers.push_back({Marker::Kind::domain, {}, -1, {}});
            }
        }
        if (at && (at->kind == "mesh" || region_index >= 0 || contact_index >= 0)) {
            for (const Pt p : named_points(out.message)) out.markers.push_back({Marker::Kind::point, {}, -1, p});
        }
    }
    return out;
}

}  // namespace NiTCAD::app::geo
