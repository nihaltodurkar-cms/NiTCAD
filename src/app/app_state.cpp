// The per-user state (app_state.hpp).
#include "NiTCAD/app/app_state.hpp"

#include <algorithm>
#include <charconv>
#include <fstream>
#include <iterator>
#include <system_error>

namespace NiTCAD::app {

namespace {

namespace fs = std::filesystem;

std::string_view trim(std::string_view s) {
    const auto b = s.find_first_not_of(" \t\r");
    if (b == std::string_view::npos) return {};
    const auto e = s.find_last_not_of(" \t\r");
    return s.substr(b, e - b + 1);
}

bool integer(std::string_view s, int& out) {
    int v = 0;
    const auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if (ec != std::errc{} || p != s.data() + s.size()) return false;
    out = v;
    return true;
}

// A path as UTF-8 and back, so names outside the system code page survive.
std::string utf8(const fs::path& p) {
    const std::u8string u = p.u8string();
    return {u.begin(), u.end()};
}
fs::path from_utf8(std::string_view s) { return fs::path(std::u8string(s.begin(), s.end())); }

}  // namespace

void AppState::add_recent(const fs::path& path) {
    remove_recent(path);
    recent.insert(recent.begin(), path);
    if (recent.size() > max_recent) recent.resize(max_recent);
}

void AppState::remove_recent(const fs::path& path) {
    std::erase_if(recent, [&](const fs::path& p) { return p == path; });
}

std::string serialize(const AppState& s) {
    std::string out = "# NiTCAD application state\n";
    for (const fs::path& p : s.recent) out += "recent = " + utf8(p) + "\n";
    const Layout& l = s.layout;
    out += "window = " + std::to_string(l.x) + " " + std::to_string(l.y) + " " + std::to_string(l.width) + " " +
           std::to_string(l.height) + "\n";
    out += "maximized = " + std::string(l.maximized ? "1" : "0") + "\n";
    out += "tree_width = " + std::to_string(l.tree_width) + "\n";
    out += "inspector_width = " + std::to_string(l.inspector_width) + "\n";
    out += "deck_height = " + std::to_string(l.deck_height) + "\n";
    return out;
}

AppState parse_state(std::string_view text) {
    AppState s;
    while (!text.empty()) {
        const auto nl = text.find('\n');
        const std::string_view line = trim(text.substr(0, nl));
        text = nl == std::string_view::npos ? std::string_view{} : text.substr(nl + 1);
        if (line.empty() || line.front() == '#') continue;
        const auto eq = line.find('=');
        if (eq == std::string_view::npos) continue;
        const std::string_view key = trim(line.substr(0, eq)), value = trim(line.substr(eq + 1));
        int v = 0;
        if (key == "recent") {
            if (!value.empty() && s.recent.size() < AppState::max_recent) {
                const fs::path p = from_utf8(value);
                if (std::ranges::find(s.recent, p) == s.recent.end()) s.recent.push_back(p);
            }
        } else if (key == "window") {
            int r[4]{};
            std::string_view rest = value;
            bool ok = true;
            for (int k = 0; k < 4 && ok; ++k) {
                rest = trim(rest);
                const auto sp = rest.find(' ');
                ok = integer(rest.substr(0, sp), r[k]) && (k == 3 || sp != std::string_view::npos);
                rest = sp == std::string_view::npos ? std::string_view{} : rest.substr(sp + 1);
            }
            if (ok && r[2] >= 0 && r[3] >= 0) s.layout = {r[0], r[1], r[2], r[3], s.layout.maximized, s.layout.tree_width,
                                                          s.layout.inspector_width, s.layout.deck_height};
        } else if (key == "maximized" && integer(value, v)) {
            s.layout.maximized = v != 0;
        } else if (key == "tree_width" && integer(value, v) && v >= 0) {
            s.layout.tree_width = v;
        } else if (key == "inspector_width" && integer(value, v) && v >= 0) {
            s.layout.inspector_width = v;
        } else if (key == "deck_height" && integer(value, v) && v >= 0) {
            s.layout.deck_height = v;
        }
    }
    return s;
}

AppState load_state(const fs::path& file) {
    std::ifstream f(file, std::ios::binary);
    if (!f) return {};
    const std::string text{std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
    return parse_state(text);
}

std::expected<void, base::Error> save_state(const AppState& state, const fs::path& file) {
    std::error_code ec;
    if (file.has_parent_path()) fs::create_directories(file.parent_path(), ec);
    std::ofstream f(file, std::ios::binary | std::ios::trunc);
    const std::string text = serialize(state);
    f.write(text.data(), static_cast<std::streamsize>(text.size()));
    if (!f) return std::unexpected(base::Error{base::ErrorCode::resource_exhausted, "cannot write " + file.string(), std::nullopt});
    return {};
}

}  // namespace NiTCAD::app
