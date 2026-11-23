// Unlock shapes (unlock.hpp; DECISIONS.md Unit 26n, DT5).
#include "NiTCAD/app/unlock.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <optional>
#include <vector>

#include "deck_text.hpp"
#include "NiTCAD/app/solid.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/mesh/generate.hpp"

namespace NiTCAD::app::unlock {

namespace {

using namespace deck_text;
using geo::Coord;
using Params = std::vector<deck::Parameter>;

constexpr std::array<std::string_view, 3> axis_names{"x", "y", "z"};

const deck::Section* last_mesh(const deck::Deck& d) {
    const deck::Section* mesh = nullptr;
    for (const auto& s : d.sections) {
        if (s.kind == "mesh") mesh = &s;
    }
    return mesh;
}

std::string mesh_type(const deck::Section& mesh) {
    const deck::Entry* t = find_entry(mesh, "type");
    return t ? std::string(trim(t->value)) : std::string("tensor");
}

// A coordinate as written and its value.
struct Word {
    std::string text;
    Coord value = 0;
};

std::optional<Word> word(std::string_view w, const Params& params) {
    w = trim(w);
    if (w.find('$') == std::string_view::npos) {
        const auto v = geo::parse_length(w);
        if (!v) return std::nullopt;
        return Word{std::string(w), *v};
    }
    const auto t = deck::substitute(w, params);
    if (!t) return std::nullopt;
    const auto v = geo::parse_length(*t);
    if (!v) return std::nullopt;
    return Word{std::string(w), *v};
}

// One tensor axis: its ends as written and its widest line spacing.
struct Axis {
    Word lo, hi;
    Coord widest = 0;
};

// The axis grammar (deck::Builder::axis): pieces joined by a top-level '+', each `linspace(a, b, n)` or a list.
std::expected<Axis, std::string> read_axis(const deck::Entry& e, const Params& params) {
    const std::string_view v = e.value;
    std::vector<std::string_view> pieces;
    std::size_t start = 0;
    int depth = 0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (v[i] == '(') ++depth;
        if (v[i] == ')') --depth;
        const bool exponent = i > 0 && (v[i - 1] == 'e' || v[i - 1] == 'E');
        if (v[i] == '+' && depth == 0 && !exponent) {
            pieces.push_back(trim(v.substr(start, i - start)));
            start = i + 1;
        }
    }
    pieces.push_back(trim(v.substr(start)));
    const auto bad = [&] { return std::unexpected("line " + std::to_string(e.line) + ": mesh lines along " + e.key + " the conversion cannot read"); };
    std::vector<Word> ends;  // each piece's first and last
    Axis a;
    for (const std::string_view piece : pieces) {
        if (piece.starts_with("linspace(") && piece.ends_with(")")) {
            const std::string_view inside = piece.substr(9, piece.size() - 10);
            std::vector<std::string_view> args;
            std::size_t from = 0;
            for (std::size_t i = 0; i <= inside.size(); ++i) {
                if (i == inside.size() || inside[i] == ',') {
                    args.push_back(trim(inside.substr(from, i - from)));
                    from = i + 1;
                }
            }
            if (args.size() != 3) return bad();
            const auto p = word(args[0], params), q = word(args[1], params);
            const auto t = deck::substitute(args[2], params);
            if (!p || !q || !t) return bad();
            char* stop = nullptr;
            const std::string count(trim(*t));
            const long n = std::strtol(count.c_str(), &stop, 10);
            if (n < 2 || stop == count.c_str() || *stop != '\0') return bad();
            a.widest = std::max(a.widest, (q->value - p->value + n - 2) / (n - 1));
            ends.push_back(*p);
            ends.push_back(*q);
        } else {
            std::optional<Word> before;
            for (const auto w : deck::value_words(piece)) {
                const auto p = word(w, params);
                if (!p) return bad();
                if (before) a.widest = std::max(a.widest, p->value - before->value);
                if (!before) ends.push_back(*p);
                before = p;
            }
            if (!before) return bad();
            ends.push_back(*before);
        }
    }
    if (ends.empty()) return bad();
    a.lo = ends.front();
    a.hi = ends.back();
    for (const Word& w : ends) {
        if (w.value < a.lo.value) a.lo = w;
        if (w.value > a.hi.value) a.hi = w;
    }
    if (a.lo.value >= a.hi.value) return bad();
    return a;
}

std::string joined(const std::vector<std::string>& w, std::string_view sep) {
    std::string out;
    for (std::size_t k = 0; k < w.size(); ++k) out += (k ? std::string(sep) : std::string()) + w[k];
    return out;
}

struct Face {
    int axis = 0;
    bool high = false;
};

std::optional<Face> face_of(std::string_view name) {
    for (int a = 0; a < 3; ++a) {
        if (name == std::string(axis_names[static_cast<std::size_t>(a)]) + "_min") return Face{a, false};
        if (name == std::string(axis_names[static_cast<std::size_t>(a)]) + "_max") return Face{a, true};
    }
    return std::nullopt;
}

// The tensor deck's text, converted; interface_resolution written when `s` is given.
struct Converted {
    std::string text;
    Coord widest = 0;
    std::vector<std::string> notes;
};

std::expected<Converted, std::string> convert_tensor(std::string_view text, std::optional<Coord> s) {
    const auto d = deck::parse(text);
    if (!d) return std::unexpected(d.error().message);
    const auto params = deck::parameters(*d);
    if (!params) return std::unexpected(params.error().message);
    const deck::Section* mesh = last_mesh(*d);
    if (!mesh || mesh_type(*mesh) != "tensor" || !find_entry(*mesh, "z")) return std::unexpected(std::string("only a 3D tensor or extrude3d deck is unlocked"));
    for (const auto& sec : d->sections) {
        if (sec.kind == "adapt") {
            return std::unexpected("line " + std::to_string(sec.line) + ": an [adapt] section refines a tensor mesh, not a layered one (remove it first)");
        }
    }
    Converted out;
    std::array<Axis, 3> ax;
    for (int k = 0; k < 3; ++k) {
        const deck::Entry* e = find_entry(*mesh, axis_names[static_cast<std::size_t>(k)]);
        if (!e) return std::unexpected(std::string("a tensor mesh needs x, then y, then z"));
        auto a = read_axis(*e, *params);
        if (!a) return std::unexpected(a.error());
        ax[static_cast<std::size_t>(k)] = std::move(*a);
        out.widest = std::max(out.widest, ax[static_cast<std::size_t>(k)].widest);
    }
    const auto lo = [&](int k) -> const Word& { return ax[static_cast<std::size_t>(k)].lo; };
    const auto hi = [&](int k) -> const Word& { return ax[static_cast<std::size_t>(k)].hi; };
    // A coordinate along axis k inside the domain: as written, else the domain's end (the box cut to it).
    bool cut = false;
    const auto inside = [&](int k, const Word& w) -> std::string {
        if (w.value < lo(k).value) return cut = true, lo(k).text;
        if (w.value > hi(k).value) return cut = true, hi(k).text;
        return w.text;
    };

    Rewrite rw(text);
    // Regions: each box kept as a box solid, cut to the domain when it reaches past it; one wholly outside dropped.
    for (const auto& sec : d->sections) {
        if (sec.kind != "region") continue;
        SectionEdit e;
        for (std::size_t k = 0; k < sec.entries.size(); ++k) {
            const auto& en = sec.entries[k];
            if (en.key != "box") continue;
            const auto w = deck::value_words(en.value);
            if (w.size() != 6) return std::unexpected("line " + std::to_string(en.line) + ": a box the conversion cannot read");
            std::vector<std::string> t;
            cut = false;
            bool empty = false;
            for (int j = 0; j < 6; ++j) {
                const auto p = word(w[static_cast<std::size_t>(j)], *params);
                if (!p) return std::unexpected("line " + std::to_string(en.line) + ": a box the conversion cannot read");
                t.push_back(inside(j / 2, *p));
            }
            for (int a = 0; a < 3; ++a) {
                const auto p = word(t[static_cast<std::size_t>(2 * a)], *params), q = word(t[static_cast<std::size_t>(2 * a + 1)], *params);
                empty = empty || !p || !q || p->value >= q->value;
            }
            if (empty) {
                e.replace[k] = std::nullopt;
                out.notes.push_back("a box of region '" + sec.name + "' outside the domain is dropped");
            } else if (cut) {
                e.replace[k] = "box = " + joined(t, " ");
                out.notes.push_back("a box of region '" + sec.name + "' is cut to the domain");
            }
        }
        if (!e.replace.empty()) rw.change(sec, rw.rewrite(sec, e));
    }
    // Contacts: the face and box become a patch box on that face, named after the contact.
    std::vector<std::string> patches;
    for (const auto& sec : d->sections) {
        if (sec.kind != "contact") continue;
        const deck::Entry* patch = find_entry(sec, "patch");
        if (!patch) return std::unexpected("[contact " + sec.name + "] has no patch");
        const auto f = face_of(trim(patch->value));
        if (!f) return std::unexpected("line " + std::to_string(patch->line) + ": a contact on '" + std::string(trim(patch->value)) + "', not a face of the domain");
        std::vector<const deck::Entry*> boxes;
        for (const auto& en : sec.entries) {
            if (en.key == "box") boxes.push_back(&en);
        }
        if (boxes.size() > 1) return std::unexpected("[contact " + sec.name + "] has " + std::to_string(boxes.size()) + " boxes: a layered deck takes one patch box a contact");
        std::array<std::string, 3> lo_t, hi_t;
        for (int a = 0; a < 3; ++a) {
            lo_t[static_cast<std::size_t>(a)] = lo(a).text;
            hi_t[static_cast<std::size_t>(a)] = hi(a).text;
        }
        if (!boxes.empty()) {
            const auto w = deck::value_words(boxes.front()->value);
            if (w.size() != 6) return std::unexpected("line " + std::to_string(boxes.front()->line) + ": a box the conversion cannot read");
            for (int a = 0; a < 3; ++a) {
                const auto p = word(w[static_cast<std::size_t>(2 * a)], *params), q = word(w[static_cast<std::size_t>(2 * a + 1)], *params);
                if (!p || !q) return std::unexpected("line " + std::to_string(boxes.front()->line) + ": a box the conversion cannot read");
                lo_t[static_cast<std::size_t>(a)] = inside(a, *p);
                hi_t[static_cast<std::size_t>(a)] = inside(a, *q);
            }
        }
        const std::string& at = f->high ? hi(f->axis).text : lo(f->axis).text;
        lo_t[static_cast<std::size_t>(f->axis)] = hi_t[static_cast<std::size_t>(f->axis)] = at;
        patches.push_back("patch " + sec.name + " = " + lo_t[0] + " " + lo_t[1] + " " + lo_t[2] + "  " + hi_t[0] + " " + hi_t[1] + " " + hi_t[2]);
        SectionEdit e;
        for (std::size_t k = 0; k < sec.entries.size(); ++k) {
            if (&sec.entries[k] == patch) e.replace[k] = "patch = " + sec.name;
            if (sec.entries[k].key == "box") e.replace[k] = std::nullopt;
        }
        rw.change(sec, rw.rewrite(sec, e));
    }
    // [mesh]: the type, layered, the domain, the patches, the sizes; the mesh lines and straddles dropped.
    {
        SectionEdit e;
        int type_at = -1;
        for (std::size_t k = 0; k < mesh->entries.size(); ++k) {
            const auto& en = mesh->entries[k];
            if (en.key == "type") {
                type_at = static_cast<int>(k);
                e.replace[k] = "type = generate3d";
            } else if (en.key == "x" || en.key == "y" || en.key == "z" || en.key.starts_with("straddle ")) {
                e.replace[k] = std::nullopt;
            }
        }
        if (type_at < 0) e.after.emplace(-1, "type = generate3d");
        e.after.emplace(type_at, "layered = on");
        e.after.emplace(type_at, "domain = " + lo(0).text + " " + hi(0).text + " " + lo(1).text + " " + hi(1).text + " " + lo(2).text + " " + hi(2).text);
        for (auto& p : patches) e.after.emplace(type_at, std::move(p));
        e.after.emplace(type_at, "max_size = " + geo::format_length(out.widest));
        if (s) e.after.emplace(type_at, "interface_resolution = " + geo::format_length(*s));
        rw.change(*mesh, rw.rewrite(*mesh, e));
    }
    out.text = rw.apply();
    return out;
}

// s (cm) rounded down to two significant figures in lattice steps (0.1 nm).
Coord tidy_resolution(double s_cm) {
    const double q = std::floor(s_cm / geo::quantum_cm);
    if (!(q >= 1.0)) return 0;
    double scale = 1.0;
    while (q / scale >= 100.0) scale *= 10.0;
    return static_cast<Coord>(std::floor(q / scale) * scale);
}

}  // namespace

Kind kind_of(std::string_view text) {
    const auto d = deck::parse(text);
    if (!d) return Kind::none;
    const deck::Section* mesh = last_mesh(*d);
    if (!mesh) return Kind::none;
    const std::string t = mesh_type(*mesh);
    if (t == "extrude3d") return Kind::extruded;
    if (t == "tensor" && find_entry(*mesh, "z")) return Kind::tensor;
    return Kind::none;
}

std::expected<Unlocked, std::string> unlock(std::string_view text) {
    switch (kind_of(text)) {
    case Kind::none: return std::unexpected(std::string("only a 3D tensor or extrude3d deck is unlocked"));
    case Kind::extruded: {
        auto m = solid::read_model(text);
        if (!m) return std::unexpected(m.error());
        if (m->mesh != solid::MeshKind::extruded) return std::unexpected(m->why);
        const auto note = solid::convert_to_layered(*m);
        if (!note) return std::unexpected(note.error());
        auto written = solid::write_model(text, *m);
        if (!written) return std::unexpected(written.error());
        return Unlocked{std::move(*written), *note, 0.0};
    }
    case Kind::tensor: break;
    }
    const auto first = convert_tensor(text, std::nullopt);
    if (!first) return std::unexpected(first.error());
    const auto m = solid::read_model(first->text);
    if (!m) return std::unexpected(m.error());
    if (m->mesh != solid::MeshKind::layered) return std::unexpected(m->why);
    const double h = static_cast<double>(first->widest) * geo::quantum_cm;
    // MT12 (1)'s rules only (no trial meshes: a 3D device meshed at every trial s takes minutes); the build meshes it.
    const auto complex = mesh::layered_partition3d(solid::to_layered(*m));
    if (!complex) return std::unexpected("the 3D generator refuses the layout: " + complex.error().message);
    mesh::MeshOptions3D o;
    o.max_size = h;
    const auto s = mesh::largest_resolution(*complex, o);
    if (!s) return std::unexpected("the 3D generator refuses the layout: " + s.error().message);
    const Coord sq = std::isfinite(*s) ? std::min(tidy_resolution(*s), first->widest) : first->widest;
    if (sq <= 0) return std::unexpected(std::format("the layout needs an interface resolution under 0.1 nm ({:.3g} nm): too fine to mesh", *s * 1e7));
    auto done = convert_tensor(text, sq);
    if (!done) return std::unexpected(done.error());
    Unlocked u;
    u.text = std::move(done->text);
    u.note = "the mesh becomes a layered generate3d mesh: the mesh lines' ends are the domain, each region box a box solid, each "
             "contact a patch box; max_size " + geo::format_length(first->widest) + ", interface_resolution " + geo::format_length(sq);
    if (const auto e = solid::extent(*m)) {
        u.nodes = solid::node_estimate(*e, h, static_cast<double>(sq) * geo::quantum_cm);
        u.note += std::format(" (about {:.0f} nodes)", u.nodes);
    }
    for (const auto& n : done->notes) u.note += "; " + n;
    return u;
}

}  // namespace NiTCAD::app::unlock
