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

namespace {

namespace mkey {
constexpr std::uint16_t format = 1, schema = 2, model_class = 3, hyper = 4, shrinkage = 5, base = 6, tree = 7,
                        material = 8, feature_min = 9, feature_max = 10, category = 11, importance = 12,
                        manifest_digest = 13, training = 14, held_out = 15, training_samples = 16, provenance = 17,
                        metrics_counts = 18, metrics = 19;
constexpr std::uint16_t t_ints = 1, t_reals = 2;      // a tree: (feature, left, right) per node; (threshold, value)
constexpr std::uint16_t e_digest = 1, e_path = 2;     // a manifest entry
constexpr std::uint16_t c_values = 1;                 // a feature's categories
}  // namespace mkey

void write_entries(Writer& w, std::uint16_t k, const std::vector<results::ManifestEntry>& entries) {
    for (const results::ManifestEntry& e : entries) {
        Writer one;
        one.u64(mkey::e_digest, e.digest);
        one.text(mkey::e_path, e.path);
        w.record(k, one);
    }
}

std::expected<std::vector<results::ManifestEntry>, base::Error> read_entries(const Reader& r, std::uint16_t k) {
    NITCAD_IO_TRY(list, r.records(k));
    std::vector<results::ManifestEntry> out;
    for (const Reader& one : list) {
        NITCAD_IO_TRY(d, one.u64(mkey::e_digest));
        NITCAD_IO_TRY(p, one.text(mkey::e_path));
        out.push_back({d, std::move(p)});
    }
    return out;
}

base::Error model_refused(const char* what, std::uint64_t found, std::uint64_t ours) {
    return invalid(std::string("a learned model of ") + what + " " + std::to_string(found) + ", this build's is " +
                   std::to_string(ours) + " (retrain the model)");
}

}  // namespace

std::expected<void, base::Error> save_learned_model(const std::filesystem::path& path, const results::LearnedModel& m) {
    namespace k = mkey;
    Writer w;
    w.u64(k::format, m.format);
    w.u64(k::schema, m.schema);
    w.text(k::model_class, m.model_class);
    w.i64s(k::hyper, std::vector<std::int64_t>{static_cast<std::int64_t>(m.tree_count), static_cast<std::int64_t>(m.depth),
                                               static_cast<std::int64_t>(m.min_leaf),
                                               static_cast<std::int64_t>(m.cut_points)});
    w.f64(k::shrinkage, m.shrinkage);
    w.f64(k::base, m.base);
    for (const auto& tree : m.trees) {
        std::vector<std::int64_t> ints;
        std::vector<double> reals;
        for (const results::TreeNode& n : tree) {
            ints.insert(ints.end(), {n.feature, n.left, n.right});
            reals.insert(reals.end(), {n.threshold, n.value});
        }
        Writer t;
        t.i64s(k::t_ints, ints);
        t.f64s(k::t_reals, reals);
        w.record(k::tree, t);
    }
    write_texts(w, k::material, m.materials);
    w.f64s(k::feature_min, m.feature_min);
    w.f64s(k::feature_max, m.feature_max);
    for (const auto& c : m.categories) {
        Writer one;
        one.f64s(k::c_values, c);
        w.record(k::category, one);
    }
    w.f64s(k::importance, m.importance);
    w.u64(k::manifest_digest, m.manifest_digest);
    write_entries(w, k::training, m.training);
    write_entries(w, k::held_out, m.held_out);
    w.u64(k::training_samples, m.training_samples);
    w.record(k::provenance, detail::write_provenance_record(m.provenance));
    w.i64s(k::metrics_counts, std::vector<std::int64_t>{static_cast<std::int64_t>(m.metrics.samples),
                                                        static_cast<std::int64_t>(m.metrics.first_pass_samples)});
    w.f64s(k::metrics, std::vector<double>{m.metrics.rank_correlation, m.metrics.ridge_rank_correlation, m.metrics.recall,
                                           m.metrics.precision});
    return detail::write_file(path, Content::learned_model, w);
}

std::expected<results::LearnedModel, base::Error> load_learned_model(const std::filesystem::path& path) {
    namespace k = mkey;
    NITCAD_IO_TRY(bytes, detail::read_file(path, Content::learned_model));
    NITCAD_IO_TRY(r, Reader::parse(bytes));
    results::LearnedModel m;
    NITCAD_IO_TRY(format, r.u64(k::format));
    if (format != results::learn_model_format) {
        return std::unexpected(model_refused("format version", format, results::learn_model_format));
    }
    NITCAD_IO_TRY(schema, r.u64(k::schema));
    if (schema != results::learn_schema_version) {
        return std::unexpected(model_refused("feature schema version", schema, results::learn_schema_version));
    }
    m.format = format;
    m.schema = schema;
    NITCAD_IO_TRY(model_class, r.text(k::model_class));
    if (model_class != "gradient_boosted_trees") return std::unexpected(corrupt("an unknown model class"));
    m.model_class = std::move(model_class);
    NITCAD_IO_TRY(hyper, r.i64s(k::hyper));
    if (hyper.size() != 4) return std::unexpected(corrupt("a learned model's hyperparameters"));
    m.tree_count = static_cast<std::uint64_t>(hyper[0]);
    m.depth = static_cast<std::uint64_t>(hyper[1]);
    m.min_leaf = static_cast<std::uint64_t>(hyper[2]);
    m.cut_points = static_cast<std::uint64_t>(hyper[3]);
    NITCAD_IO_TRY(shrinkage, r.f64(k::shrinkage));
    m.shrinkage = shrinkage;
    NITCAD_IO_TRY(base, r.f64(k::base));
    m.base = base;
    NITCAD_IO_TRY(trees, r.records(k::tree));
    for (const Reader& tr : trees) {
        NITCAD_IO_TRY(ints, tr.i64s(k::t_ints));
        NITCAD_IO_TRY(reals, tr.f64s(k::t_reals));
        if (ints.empty() || ints.size() % 3 != 0 || reals.size() != ints.size() / 3 * 2) {
            return std::unexpected(corrupt("a learned model's tree"));
        }
        const std::size_t nodes = ints.size() / 3;
        std::vector<results::TreeNode> tree(nodes);
        for (std::size_t i = 0; i < nodes; ++i) {
            results::TreeNode& n = tree[i];
            n.feature = ints[3 * i];
            n.left = ints[3 * i + 1];
            n.right = ints[3 * i + 2];
            n.threshold = reals[2 * i];
            n.value = reals[2 * i + 1];
            // A child always follows its parent, so traversal ends.
            const auto ok = [&](std::int64_t c) { return c > static_cast<std::int64_t>(i) && c < static_cast<std::int64_t>(nodes); };
            if (n.feature >= 0 && (n.feature >= static_cast<std::int64_t>(F) || !ok(n.left) || !ok(n.right))) {
                return std::unexpected(corrupt("a learned model's tree node"));
            }
            if (n.feature < -1) return std::unexpected(corrupt("a learned model's tree node"));
        }
        m.trees.push_back(std::move(tree));
    }
    if (m.trees.size() != m.tree_count) return std::unexpected(corrupt("a learned model's tree count"));
    NITCAD_IO_TRY(materials, read_texts(r, k::material));
    m.materials = std::move(materials);
    NITCAD_IO_TRY(fmin, r.f64s(k::feature_min));
    NITCAD_IO_TRY(fmax, r.f64s(k::feature_max));
    NITCAD_IO_TRY(importance, r.f64s(k::importance));
    if (fmin.size() != F || fmax.size() != F || importance.size() != F) {
        return std::unexpected(corrupt("a learned model's per-feature arrays"));
    }
    m.feature_min = std::move(fmin);
    m.feature_max = std::move(fmax);
    m.importance = std::move(importance);
    NITCAD_IO_TRY(categories, r.records(k::category));
    if (categories.size() != F) return std::unexpected(corrupt("a learned model's categories"));
    for (const Reader& c : categories) {
        NITCAD_IO_TRY(values, c.f64s_or_empty(k::c_values));
        m.categories.push_back(std::move(values));
    }
    NITCAD_IO_TRY(md, r.u64(k::manifest_digest));
    m.manifest_digest = md;
    NITCAD_IO_TRY(training, read_entries(r, k::training));
    m.training = std::move(training);
    NITCAD_IO_TRY(held_out, read_entries(r, k::held_out));
    m.held_out = std::move(held_out);
    NITCAD_IO_TRY(samples, r.u64(k::training_samples));
    m.training_samples = samples;
    NITCAD_IO_TRY(prov, r.record(k::provenance));
    NITCAD_IO_TRY(provenance, detail::read_provenance_record(prov));
    m.provenance = std::move(provenance);
    NITCAD_IO_TRY(counts, r.i64s(k::metrics_counts));
    NITCAD_IO_TRY(metrics, r.f64s(k::metrics));
    if (counts.size() != 2 || metrics.size() != 4) return std::unexpected(corrupt("a learned model's metrics"));
    m.metrics = {static_cast<std::uint64_t>(counts[0]), static_cast<std::uint64_t>(counts[1]), metrics[0], metrics[1],
                 metrics[2], metrics[3]};
    return m;
}

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
