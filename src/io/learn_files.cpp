// Cross-run learning files (Unit 17b; io/learn_files.hpp). Keys are local to each content kind.
#include "NiTCAD/io/learn_files.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "container.hpp"
#include "record.hpp"

namespace NiTCAD::io {

namespace {

using detail::corrupt;
using detail::Reader;
using detail::Writer;

constexpr std::size_t F = results::learn_feature_count;
constexpr const char* manifest_magic = "nitcad-learn-manifest";
constexpr int manifest_version = 1;

namespace hkey {
constexpr std::uint16_t schema = 1, deck = 2, run = 3, material = 4, pass = 5;
constexpr std::uint16_t p_ints = 1, p_features = 2, p_labels = 3;
}  // namespace hkey

base::Error invalid(std::string message) { return detail::invalid(std::move(message)); }

base::Error schema_refused(std::uint64_t found) {
    return invalid("learning schema version " + std::to_string(found) + ", this build's is " +
                   std::to_string(results::learn_schema_version) + " (retrain from new histories)");
}

std::string hex16(std::uint64_t v) {
    char b[17];
    std::snprintf(b, sizeof b, "%016llx", static_cast<unsigned long long>(v));
    return b;
}

void write_texts(Writer& w, std::uint16_t k, const std::vector<std::string>& texts) {
    for (const std::string& t : texts) {
        Writer one;
        one.text(1, t);
        w.record(k, one);
    }
}

std::expected<std::vector<std::string>, base::Error> read_texts(const Reader& r, std::uint16_t k) {
    NITCAD_IO_TRY(list, r.records(k));
    std::vector<std::string> out;
    for (const Reader& one : list) {
        NITCAD_IO_TRY(t, one.text(1));
        out.push_back(std::move(t));
    }
    return out;
}

}  // namespace

std::expected<void, base::Error> save_adapt_history(const std::filesystem::path& path, const results::AdaptHistory& h) {
    namespace k = hkey;
    if (h.schema != results::learn_schema_version) return std::unexpected(schema_refused(h.schema));
    Writer w;
    w.u64(k::schema, h.schema);
    w.text(k::deck, h.deck_text);
    w.record(k::run, detail::write_run_record(h.run));
    write_texts(w, k::material, h.materials);
    for (const results::HistoryPass& p : h.passes) {
        if (p.features.size() != p.rows * F || p.labels.size() != p.rows) {
            return std::unexpected(invalid("a history pass needs rows x features and one label per row"));
        }
        Writer pw;
        pw.i64s(k::p_ints, std::vector<std::int64_t>{static_cast<std::int64_t>(p.device_digest),
                                                     static_cast<std::int64_t>(p.mesh_digest),
                                                     static_cast<std::int64_t>(p.rows)});
        pw.f64s(k::p_features, p.features);
        pw.f64s(k::p_labels, p.labels);
        w.record(k::pass, pw);
    }
    return detail::write_file(path, Content::adapt_history, w);
}

std::expected<results::AdaptHistory, base::Error> load_adapt_history(const std::filesystem::path& path) {
    namespace k = hkey;
    NITCAD_IO_TRY(bytes, detail::read_file(path, Content::adapt_history));
    NITCAD_IO_TRY(r, Reader::parse(bytes));
    results::AdaptHistory h;
    NITCAD_IO_TRY(schema, r.u64(k::schema));
    if (schema != results::learn_schema_version) return std::unexpected(schema_refused(schema));
    h.schema = schema;
    NITCAD_IO_TRY(deck, r.text(k::deck));
    h.deck_text = std::move(deck);
    NITCAD_IO_TRY(run_record, r.record(k::run));
    NITCAD_IO_TRY(run, detail::read_run_record(run_record));
    h.run = std::move(run);
    NITCAD_IO_TRY(materials, read_texts(r, k::material));
    h.materials = std::move(materials);
    NITCAD_IO_TRY(passes, r.records(k::pass));
    for (const Reader& pr : passes) {
        results::HistoryPass p;
        NITCAD_IO_TRY(ints, pr.i64s(k::p_ints));
        if (ints.size() != 3) return std::unexpected(corrupt("a history pass's header"));
        p.device_digest = static_cast<std::uint64_t>(ints[0]);
        p.mesh_digest = static_cast<std::uint64_t>(ints[1]);
        p.rows = static_cast<std::uint64_t>(ints[2]);
        NITCAD_IO_TRY(features, pr.f64s(k::p_features));
        NITCAD_IO_TRY(labels, pr.f64s(k::p_labels));
        if (features.size() != p.rows * F || labels.size() != p.rows) {
            return std::unexpected(corrupt("a history pass's features or labels"));
        }
        for (std::size_t e = 0; e < p.rows; ++e) {
            const double m = features[e * F + 1];
            if (m >= 0.0 && !(m < static_cast<double>(h.materials.size()))) {
                return std::unexpected(corrupt("a history sample's material index"));
            }
        }
        p.features = std::move(features);
        p.labels = std::move(labels);
        h.passes.push_back(std::move(p));
    }
    return h;
}

std::expected<std::uint64_t, base::Error> file_digest(const std::filesystem::path& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return std::unexpected(invalid("cannot open '" + path.string() + "'"));
    std::vector<char> raw((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    return detail::fnv1a(std::as_bytes(std::span<const char>(raw)));
}

std::expected<results::Manifest, base::Error> make_manifest(std::span<const std::filesystem::path> histories) {
    results::Manifest m;
    for (const std::filesystem::path& p : histories) {
        if (auto h = load_adapt_history(p); !h) {
            base::Error e = h.error();
            e.message = "'" + p.string() + "': " + e.message;
            return std::unexpected(std::move(e));
        }
        NITCAD_IO_TRY(d, file_digest(p));
        m.entries.push_back({d, p.generic_string()});
    }
    std::ranges::stable_sort(m.entries, {}, &results::ManifestEntry::digest);
    for (std::size_t i = 1; i < m.entries.size(); ++i) {
        if (m.entries[i].digest == m.entries[i - 1].digest) {
            return std::unexpected(invalid("'" + m.entries[i].path + "' duplicates '" + m.entries[i - 1].path + "'"));
        }
    }
    return m;
}

std::expected<void, base::Error> save_manifest(const std::filesystem::path& path, const results::Manifest& m) {
    std::string text = std::string(manifest_magic) + " " + std::to_string(manifest_version) + "\nschema " +
                       std::to_string(m.schema) + "\n";
    for (const results::ManifestEntry& e : m.entries) {
        if (e.path.find('\n') != std::string::npos) return std::unexpected(invalid("a manifest path holds a line break"));
        text += hex16(e.digest) + " " + e.path + "\n";
    }
    std::filesystem::path temporary = path;
    temporary += ".tmp";
    {
        std::ofstream f(temporary, std::ios::binary | std::ios::trunc);
        if (!f) return std::unexpected(base::Error{base::ErrorCode::resource_exhausted, "cannot write '" + path.string() + "'", std::nullopt});
        f.write(text.data(), static_cast<std::streamsize>(text.size()));
        if (!f) return std::unexpected(base::Error{base::ErrorCode::resource_exhausted, "cannot write '" + path.string() + "'", std::nullopt});
    }
    std::error_code ec;
    std::filesystem::rename(temporary, path, ec);
    if (ec) return std::unexpected(base::Error{base::ErrorCode::resource_exhausted, "cannot write '" + path.string() + "': " + ec.message(), std::nullopt});
    return {};
}

std::expected<results::Manifest, base::Error> load_manifest(const std::filesystem::path& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return std::unexpected(invalid("cannot open '" + path.string() + "'"));
    const auto bad = [&](std::size_t line, const std::string& what) {
        return invalid("'" + path.string() + "' line " + std::to_string(line) + ": " + what);
    };
    std::string line;
    if (!std::getline(f, line) || line != std::string(manifest_magic) + " " + std::to_string(manifest_version)) {
        return std::unexpected(bad(1, "not a version " + std::to_string(manifest_version) + " learning manifest"));
    }
    results::Manifest m;
    if (!std::getline(f, line) || line.rfind("schema ", 0) != 0) return std::unexpected(bad(2, "expected 'schema <n>'"));
    try {
        m.schema = std::stoull(line.substr(7));
    } catch (...) {
        return std::unexpected(bad(2, "expected 'schema <n>'"));
    }
    if (m.schema != results::learn_schema_version) return std::unexpected(schema_refused(m.schema));
    for (std::size_t n = 3; std::getline(f, line); ++n) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line.size() < 18 || line[16] != ' ') return std::unexpected(bad(n, "expected '<digest> <path>'"));
        results::ManifestEntry e;
        try {
            std::size_t used = 0;
            e.digest = std::stoull(line.substr(0, 16), &used, 16);
            if (used != 16) return std::unexpected(bad(n, "a digest is 16 hexadecimal digits"));
        } catch (...) {
            return std::unexpected(bad(n, "a digest is 16 hexadecimal digits"));
        }
        e.path = line.substr(17);
        if (!m.entries.empty() && !(m.entries.back().digest < e.digest)) {
            return std::unexpected(bad(n, "entries must be sorted by digest, without duplicates"));
        }
        m.entries.push_back(std::move(e));
    }
    return m;
}

std::expected<results::TrainingSet, base::Error> build_training_set(const results::Manifest& manifest,
                                                                    const std::filesystem::path& base,
                                                                    std::size_t max_per_history) {
    if (manifest.schema != results::learn_schema_version) return std::unexpected(schema_refused(manifest.schema));
    if (max_per_history == 0) return std::unexpected(invalid("the per-history cap must be at least 1"));
    std::vector<results::AdaptHistory> histories;
    results::TrainingSet t;
    t.devices = manifest.entries;
    for (const results::ManifestEntry& e : manifest.entries) {
        std::filesystem::path p(e.path);
        if (p.is_relative() && !base.empty()) p = base / p;
        NITCAD_IO_TRY(d, file_digest(p));
        if (d != e.digest) {
            return std::unexpected(invalid("'" + e.path + "' has digest " + hex16(d) + ", the manifest says " + hex16(e.digest)));
        }
        NITCAD_IO_TRY(h, load_adapt_history(p));
        t.materials.insert(t.materials.end(), h.materials.begin(), h.materials.end());
        histories.push_back(std::move(h));
    }
    std::ranges::sort(t.materials);
    t.materials.erase(std::ranges::unique(t.materials).begin(), t.materials.end());
    for (std::size_t dv = 0; dv < histories.size(); ++dv) {
        const results::AdaptHistory& h = histories[dv];
        std::vector<double> remap;
        for (const std::string& key : h.materials) {
            remap.push_back(static_cast<double>(std::ranges::lower_bound(t.materials, key) - t.materials.begin()));
        }
        std::size_t labelled = 0;
        for (const results::HistoryPass& p : h.passes) {
            for (const double l : p.labels) labelled += std::isnan(l) ? 0 : 1;
        }
        const std::size_t stride = (labelled + max_per_history - 1) / max_per_history;
        std::size_t seen = 0;
        for (std::size_t pi = 0; pi < h.passes.size(); ++pi) {
            const results::HistoryPass& p = h.passes[pi];
            for (std::size_t e = 0; e < p.rows; ++e) {
                if (std::isnan(p.labels[e])) continue;
                if (seen++ % std::max<std::size_t>(stride, 1) != 0) continue;
                const double* f = &p.features[e * F];
                for (std::size_t k = 0; k < F; ++k) {
                    double v = f[k];
                    if (k == 1 && v >= 0.0) v = remap[static_cast<std::size_t>(v)];
                    t.features.push_back(v);
                }
                t.labels.push_back(p.labels[e]);
                t.device.push_back(static_cast<std::uint32_t>(dv));
                t.pass.push_back(static_cast<std::uint32_t>(pi));
                ++t.rows;
            }
        }
    }
    return t;
}

}  // namespace NiTCAD::io
