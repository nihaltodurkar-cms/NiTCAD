#include "NiTCAD/app/new_device.hpp"

#include <algorithm>
#include <fstream>
#include <limits>
#include <system_error>

#include "NiTCAD/app/dimensions.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"

namespace NiTCAD::app::devices {

namespace fs = std::filesystem;

namespace {

std::string read_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return {};
    return s.substr(b, s.find_last_not_of(" \t\r\n") - b + 1);
}

// The deck with each [output] file named <stem>.<its extension>, each line rewritten in place (its key and comment
// kept), so two decks made from one family in one folder write their own results; the text as it is when it does
// not parse.
std::string with_outputs_named(std::string_view deck_text, const std::string& stem) {
    const auto d = deck::parse(deck_text);
    if (!d) return std::string(deck_text);
    std::vector<std::string> lines;
    for (std::size_t at = 0;;) {
        const auto nl = deck_text.find('\n', at);
        lines.emplace_back(deck_text.substr(at, nl == std::string_view::npos ? std::string_view::npos : nl - at));
        if (nl == std::string_view::npos) break;
        at = nl + 1;
    }
    for (const deck::Section& s : d->sections) {
        if (s.kind != "output") continue;
        for (const deck::Entry& e : s.entries) {
            const auto k = static_cast<std::size_t>(e.line - 1);
            if (k < lines.size()) lines[k] = dims::with_value(lines[k], stem + fs::path(e.value).extension().string());
        }
    }
    std::string out;
    for (std::size_t k = 0; k < lines.size(); ++k) out += (k ? "\n" : "") + lines[k];
    return out;
}

}  // namespace

std::vector<Family> list_families(const fs::path& dir) {
    std::vector<Family> out;
    for (Example& e : list_examples(dir)) {
        int dimension = 0;
        if (const auto d = deck::parse(read_file(e.path))) {
            if (const auto plan = deck::build(*d)) dimension = plan->device.mesh().dimension();
        }
        out.push_back({std::move(e), dimension});
    }
    return out;
}

std::vector<std::string> main_dimensions(std::string_view deck_text) {
    const auto rows = dims::read_dimensions(deck_text);
    std::vector<std::string> out;
    if (!rows || rows->empty()) return out;
    const std::string first = rows->front().group;
    for (const dims::Dimension& r : *rows) {
        if (r.group != first) break;
        out.push_back(r.name);
    }
    return out;
}

Thumbnail thumbnail(std::string_view deck_text) {
    Thumbnail t;
    const auto d = deck::parse(deck_text);
    if (!d) {
        t.error = d.error().message;
        return t;
    }
    const auto plan = deck::build(*d);
    if (!plan) {
        t.error = plan.error().message;
        return t;
    }
    const auto s = deck::substituted(*d);
    if (!s) {
        t.error = s.error().message;
        return t;
    }
    const auto& mesh = plan->device.mesh();
    const int dim = mesh.dimension();
    double lo[3], hi[3];
    for (int k = 0; k < 3; ++k) {
        lo[k] = std::numeric_limits<double>::max();
        hi[k] = std::numeric_limits<double>::lowest();
    }
    for (const auto& p : mesh.points()) {
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], p[static_cast<std::size_t>(k)]);
            hi[k] = std::max(hi[k], p[static_cast<std::size_t>(k)]);
        }
    }
    t.x0 = lo[0];
    t.x1 = hi[0];
    t.y0 = dim == 1 ? 0.0 : lo[1];
    t.y1 = dim == 1 ? 1.0 : hi[1];
    const double zm = 0.5 * (lo[2] + hi[2]);
    for (const deck::Section& sec : s->sections) {
        if (sec.kind != "region") continue;
        std::string material;
        std::vector<std::vector<double>> boxes;
        for (const deck::Entry& e : sec.entries) {
            if (e.key == "material") material = e.value;
            if (e.key == "box") boxes.push_back(dims::lengths(e.value));
        }
        if (boxes.empty()) {
            t.boxes.push_back({t.x0, t.y0, t.x1, t.y1, material});
            continue;
        }
        for (const auto& b : boxes) {
            if (dim == 1 && b.size() == 2) t.boxes.push_back({b[0], t.y0, b[1], t.y1, material});
            if (dim == 2 && b.size() == 4) t.boxes.push_back({b[0], b[2], b[1], b[3], material});
            if (dim == 3 && b.size() == 6 && b[4] <= zm && zm <= b[5]) t.boxes.push_back({b[0], b[2], b[1], b[3], material});
        }
    }
    return t;
}

std::expected<fs::path, std::string> create_device(std::string_view deck_text, const fs::path& folder, std::string name) {
    name = trim(std::move(name));
    if (name.empty()) return std::unexpected("Type a name for the deck.");
    if (name.find_first_of(R"(\/:*?"<>|)") != std::string::npos) return std::unexpected(R"(A name cannot hold \ / : * ? " < > |.)");
    if (!name.ends_with(".deck")) name += ".deck";
    std::error_code ec;
    fs::create_directories(folder, ec);
    if (ec) return std::unexpected("The folder " + folder.string() + " cannot be made: " + ec.message());
    const fs::path to = folder / name;
    if (fs::exists(to, ec)) return std::unexpected(name + " is already in " + folder.string() + ": choose another name");
    const std::string text = with_outputs_named(deck_text, fs::path(name).stem().string());
    std::ofstream f(to, std::ios::binary);
    f.write(text.data(), static_cast<std::streamsize>(text.size()));
    f.close();
    if (!f) return std::unexpected("The deck cannot be written to " + to.string() + ".");
    return to;
}

}  // namespace NiTCAD::app::devices
