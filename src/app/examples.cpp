// The example decks (examples.hpp).
#include "NiTCAD/app/examples.hpp"

#include <algorithm>
#include <fstream>
#include <system_error>

namespace NiTCAD::app {

namespace {

namespace fs = std::filesystem;

// A line's comment text: after '#' and blanks, trimmed; empty when the line is not a comment.
std::string comment_text(std::string line) {
    while (!line.empty() && (line.back() == '\r' || line.back() == ' ')) line.pop_back();
    const auto b = line.find_first_not_of(' ');
    if (b == std::string::npos || line[b] != '#') return {};
    const auto t = line.find_first_not_of(" #", b);
    return t == std::string::npos ? std::string() : line.substr(t);
}

bool has_decks(const fs::path& dir) {
    std::error_code ec;
    if (!fs::is_directory(dir, ec)) return false;
    for (const auto& e : fs::directory_iterator(dir, ec)) {
        if (e.path().extension() == ".deck") return true;
    }
    return false;
}

}  // namespace

std::vector<Example> list_examples(const fs::path& directory) {
    std::vector<Example> out;
    std::error_code ec;
    if (!fs::is_directory(directory, ec)) return out;
    for (const auto& e : fs::directory_iterator(directory, ec)) {
        if (!e.is_regular_file(ec) || e.path().extension() != ".deck") continue;
        Example x;
        x.path = e.path();
        std::ifstream f(e.path(), std::ios::binary);
        std::string first, second;
        std::getline(f, first);
        std::getline(f, second);
        x.title = comment_text(first);
        if (x.title.empty()) x.title = e.path().stem().string();
        else x.description = comment_text(second);
        out.push_back(std::move(x));
    }
    std::ranges::sort(out, {}, [](const Example& x) { return x.path.filename(); });
    return out;
}

fs::path examples_directory(const fs::path& executable_dir, const fs::path& fallback) {
    const fs::path beside = executable_dir / "examples";
    return has_decks(beside) ? beside : fallback;
}

std::expected<fs::path, base::Error> copy_example(const Example& example, const fs::path& to) {
    std::error_code ec;
    fs::create_directories(to, ec);
    if (ec) {
        return std::unexpected(base::Error{base::ErrorCode::resource_exhausted,
                                           "cannot create " + to.string() + ": " + ec.message(), std::nullopt});
    }
    const fs::path copy = to / example.path.filename();
    if (fs::exists(copy, ec)) return copy;
    fs::copy_file(example.path, copy, ec);
    if (ec) {
        return std::unexpected(base::Error{base::ErrorCode::resource_exhausted,
                                           "cannot copy " + example.path.string() + ": " + ec.message(), std::nullopt});
    }
    return copy;
}

}  // namespace NiTCAD::app
