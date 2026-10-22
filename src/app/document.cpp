// The deck document (document.hpp).
#include "NiTCAD/app/document.hpp"

#include <charconv>
#include <fstream>
#include <iterator>
#include <system_error>
#include <utility>

#include "NiTCAD/deck/syntax.hpp"

namespace NiTCAD::app {

namespace {

namespace fs = std::filesystem;

std::string source_of(const fs::path& path) { return path.empty() ? std::string("<untitled>") : path.filename().string(); }

base::Error write_error(const fs::path& path, std::string_view what) {
    return {base::ErrorCode::resource_exhausted, "cannot write " + path.string() + ": " + std::string(what), std::nullopt};
}

}  // namespace

int error_line(const base::Error& error, std::string_view source) {
    const std::string_view m = error.message;
    if (!m.starts_with(source) || m.size() <= source.size() || m[source.size()] != ':') return 0;
    const std::string_view rest = m.substr(source.size() + 1);
    int line = 0;
    const auto [end, ec] = std::from_chars(rest.data(), rest.data() + rest.size(), line);
    if (ec != std::errc{} || end == rest.data() || end == rest.data() + rest.size() || *end != ':') return 0;
    return line > 0 ? line : 0;
}

CheckResult check_deck(std::string_view text, const fs::path& path) {
    CheckResult out;
    const std::string source = source_of(path);
    auto parsed = deck::parse(text, source, path.empty() ? fs::path{} : path.parent_path());
    if (!parsed) {
        out.diagnostics.push_back({error_line(parsed.error(), source), parsed.error().message});
        return out;
    }
    auto plan = deck::build(*parsed);
    if (!plan) {
        out.diagnostics.push_back({error_line(plan.error(), source), plan.error().message});
        return out;
    }
    out.plan = std::move(*plan);
    return out;
}

DeckDocument DeckDocument::untitled(std::string text) {
    DeckDocument d;
    d.text_ = std::move(text);
    return d;
}

std::expected<DeckDocument, base::Error> DeckDocument::open(const fs::path& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return std::unexpected(base::Error{base::ErrorCode::invalid_input, "cannot read " + path.string(), std::nullopt});
    DeckDocument d;
    d.text_.assign(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
    if (f.bad()) return std::unexpected(base::Error{base::ErrorCode::invalid_input, "cannot read " + path.string(), std::nullopt});
    d.path_ = path;
    return d;
}

std::string DeckDocument::title() const {
    std::string t = path_.empty() ? std::string("Untitled") : path_.filename().string();
    if (modified_) t += " *";
    return t;
}

void DeckDocument::set_text(std::string text) {
    if (text == text_) return;
    text_ = std::move(text);
    modified_ = true;
}

std::expected<void, base::Error> DeckDocument::save() {
    if (path_.empty()) {
        return std::unexpected(base::Error{base::ErrorCode::invalid_input, "an untitled deck needs a file name", std::nullopt});
    }
    // A temporary file next to the deck, then renamed over it: the old file stays whole if the write
    // fails.
    fs::path temporary = path_;
    temporary += ".saving";
    {
        std::ofstream f(temporary, std::ios::binary | std::ios::trunc);
        if (!f) return std::unexpected(write_error(path_, "cannot create the file"));
        f.write(text_.data(), static_cast<std::streamsize>(text_.size()));
        f.flush();
        if (!f) {
            f.close();
            std::error_code ignored;
            fs::remove(temporary, ignored);
            return std::unexpected(write_error(path_, "the write failed"));
        }
    }
    std::error_code ec;
    fs::rename(temporary, path_, ec);
    if (ec) {
        std::error_code ignored;
        fs::remove(temporary, ignored);
        return std::unexpected(write_error(path_, ec.message()));
    }
    modified_ = false;
    return {};
}

std::expected<void, base::Error> DeckDocument::save_as(const fs::path& path) {
    const fs::path old = std::exchange(path_, path);
    if (auto r = save(); !r) {
        path_ = old;
        return r;
    }
    return {};
}

}  // namespace NiTCAD::app
