// `nitcad learn` (Unit 17b; learn.hpp).
#include "learn.hpp"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <numeric>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include "NiTCAD/io/learn_files.hpp"
#include "NiTCAD/results/learn.hpp"
#include "NiTCAD/solve/adapt_learn.hpp"

namespace NiTCAD::cli {
namespace {

int usage() {
    std::fputs("usage: nitcad learn manifest <manifest> <history.nithist>...\n"
               "       nitcad learn train <manifest> <model.nitlearn> [--holdout <history.nithist>]...\n",
               stderr);
    return 64;
}

int fail(const base::Error& e) {
    std::fprintf(stderr, "nitcad learn: %s\n", e.message.c_str());
    return 1;
}

int manifest(int argc, char** argv) {
    if (argc < 5) return usage();  // nitcad learn manifest <manifest> <history>...
    std::vector<std::filesystem::path> histories;
    for (int i = 4; i < argc; ++i) histories.emplace_back(argv[i]);
    auto m = io::make_manifest(histories);
    if (!m) return fail(m.error());
    if (auto w = io::save_manifest(argv[3], *m); !w) return fail(w.error());
    for (const results::ManifestEntry& e : m->entries) {
        std::printf("%016llx %s\n", static_cast<unsigned long long>(e.digest), e.path.c_str());
    }
    std::printf("manifest: %zu histories, schema %llu\n", m->entries.size(),
                static_cast<unsigned long long>(m->schema));
    return 0;
}

// The manifest entry a --holdout path names: the same text, or the same file.
std::ptrdiff_t entry_of(const results::Manifest& m, const std::filesystem::path& base, const std::string& arg) {
    for (std::size_t i = 0; i < m.entries.size(); ++i) {
        if (m.entries[i].path == arg) return static_cast<std::ptrdiff_t>(i);
        std::filesystem::path p(m.entries[i].path);
        if (p.is_relative()) p = base / p;
        std::error_code ec;
        if (std::filesystem::equivalent(p, arg, ec) && !ec) return static_cast<std::ptrdiff_t>(i);
    }
    return -1;
}

int train(int argc, char** argv) {
    if (argc < 5) return usage();  // nitcad learn train <manifest> <model> [--holdout <history>]...
    const std::filesystem::path manifest_path = argv[3], model_path = argv[4];
    auto m = io::load_manifest(manifest_path);
    if (!m) return fail(m.error());
    const std::filesystem::path base = manifest_path.parent_path();
    std::vector<std::size_t> held;
    for (int i = 5; i < argc; ++i) {
        if (std::string_view(argv[i]) != "--holdout" || i + 1 >= argc) return usage();
        const std::ptrdiff_t e = entry_of(*m, base, argv[++i]);
        if (e < 0) {
            return fail({base::ErrorCode::invalid_input, std::string("'") + argv[i] + "' is not in the manifest", std::nullopt});
        }
        held.push_back(static_cast<std::size_t>(e));
    }
    if (held.empty()) held = solve::default_held_out(m->entries.size());
    auto set = io::build_training_set(*m, base);
    if (!set) return fail(set.error());
    auto model = solve::train_predictor(*set, held);
    if (!model) return fail(model.error());
    if (auto w = io::save_learned_model(model_path, *model); !w) return fail(w.error());
    const results::LearnedModel& r = *model;
    std::printf("trained: %zu trees of depth %llu on %llu samples from %zu devices; held out %zu devices\n",
                r.trees.size(), static_cast<unsigned long long>(r.depth),
                static_cast<unsigned long long>(r.training_samples), r.training.size(), r.held_out.size());
    for (const results::ManifestEntry& e : r.held_out) std::printf("held out: %s\n", e.path.c_str());
    std::printf("held-out rank correlation: trees %.4f, ridge %.4f (%llu samples)\n", r.metrics.rank_correlation,
                r.metrics.ridge_rank_correlation, static_cast<unsigned long long>(r.metrics.samples));
    std::printf("held-out first pass: recall %.4f, precision %.4f (%llu samples)\n", r.metrics.recall,
                r.metrics.precision, static_cast<unsigned long long>(r.metrics.first_pass_samples));
    std::vector<std::size_t> order(r.importance.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::ranges::stable_sort(order, [&](std::size_t a, std::size_t b) { return r.importance[a] > r.importance[b]; });
    double total = 0.0;
    for (const double v : r.importance) total += v;
    for (const std::size_t k : order) {
        if (!(r.importance[k] > 0.0)) break;
        std::printf("importance: %-24s %.4f\n", std::string(results::learn_feature_names[k]).c_str(),
                    total > 0.0 ? r.importance[k] / total : 0.0);
    }
    return 0;
}

}  // namespace

int learn(int argc, char** argv) {
    if (argc < 3) return usage();
    const std::string_view sub = argv[2];
    if (sub == "manifest") return manifest(argc, argv);
    if (sub == "train") return train(argc, argv);
    return usage();
}

}  // namespace NiTCAD::cli
