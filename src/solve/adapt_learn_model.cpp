// Cross-run learning: the model (Unit 17b, slice 2; DECISIONS.md LB5, LB7, LB10; ID9, ID10).
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <string>
#include <utility>

#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/solve/adapt_learn.hpp"
#include "NiTCAD/solve/provenance.hpp"

namespace NiTCAD::solve {
namespace {

base::Error invalid(std::string message) { return {base::ErrorCode::invalid_input, std::move(message), std::nullopt}; }

constexpr std::size_t F = results::learn_feature_count;

// FNV-1a 64 over the entries' digests (the data, not where it lies).
std::uint64_t manifest_digest_of(std::span<const results::ManifestEntry> entries) {
    std::uint64_t h = 0xcbf29ce484222325ull;
    for (const results::ManifestEntry& e : entries) {
        for (int i = 0; i < 8; ++i) {
            h ^= (e.digest >> (8 * i)) & 0xffu;
            h *= 0x100000001b3ull;
        }
    }
    return h;
}

// Ranks, averaged over ties (1-based).
std::vector<double> ranks_of(std::span<const double> v) {
    std::vector<std::size_t> order(v.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::ranges::stable_sort(order, [&](std::size_t a, std::size_t b) { return v[a] < v[b]; });
    std::vector<double> r(v.size());
    for (std::size_t i = 0; i < order.size();) {
        std::size_t j = i;
        while (j + 1 < order.size() && v[order[j + 1]] == v[order[i]]) ++j;
        const double avg = 0.5 * static_cast<double>(i + j) + 1.0;
        for (std::size_t k = i; k <= j; ++k) r[order[k]] = avg;
        i = j + 1;
    }
    return r;
}

std::size_t leaf_of(const std::vector<results::TreeNode>& tree, const double* row) {
    std::size_t node = 0;
    while (tree[node].feature >= 0) {
        const results::TreeNode& t = tree[node];
        node = static_cast<std::size_t>(row[t.feature] <= t.threshold ? t.left : t.right);
    }
    return node;
}

// The ridge's design row (intercept first) for a feature row.
void ridge_row(const Ridge& r, const double* row, std::vector<double>& x) {
    x.clear();
    x.push_back(1.0);
    for (std::size_t k = 0; k < F; ++k) {
        if (results::learn_feature_is_category(k)) {
            for (const double c : r.categories[k]) x.push_back(row[k] == c ? 1.0 : 0.0);
        } else if (r.scale[k] > 0.0) {
            x.push_back((row[k] - r.mean[k]) / r.scale[k]);
        }
    }
}

}  // namespace

std::vector<double> cut_points(std::span<const double> values, std::size_t max_cuts) {
    std::vector<double> v(values.begin(), values.end());
    std::ranges::sort(v);
    std::vector<double> distinct = v;
    distinct.erase(std::ranges::unique(distinct).begin(), distinct.end());
    std::vector<double> cuts;
    if (distinct.size() <= 1 || max_cuts == 0) return cuts;
    if (distinct.size() <= max_cuts) {
        cuts.assign(distinct.begin(), distinct.end() - 1);
        return cuts;
    }
    const std::size_t n = v.size();
    for (std::size_t j = 1; j <= max_cuts; ++j) {
        const double c = v[std::min(n - 1, (j * n) / (max_cuts + 1))];
        if (c == distinct.back()) continue;
        if (cuts.empty() || c != cuts.back()) cuts.push_back(c);
    }
    return cuts;
}

BoostedTrees fit_boosted_trees(std::span<const double> features, std::span<const double> labels,
                               std::span<const std::size_t> rows, const TreeOptions& o) {
    BoostedTrees out;
    out.importance.assign(F, 0.0);
    const std::size_t n = rows.size();
    if (n == 0) return out;
    // Cut points and each row's bin per feature (bin b: x <= cuts[b]; the last bin above every cut).
    std::vector<std::vector<double>> cuts(F);
    std::vector<std::vector<std::uint8_t>> bin(F, std::vector<std::uint8_t>(n));
    {
        std::vector<double> column(n);
        for (std::size_t k = 0; k < F; ++k) {
            for (std::size_t i = 0; i < n; ++i) column[i] = features[rows[i] * F + k];
            cuts[k] = cut_points(column, std::min<std::size_t>(o.cut_points, 254));
            for (std::size_t i = 0; i < n; ++i) {
                bin[k][i] = static_cast<std::uint8_t>(std::ranges::lower_bound(cuts[k], column[i]) - cuts[k].begin());
            }
        }
    }
    double sum = 0.0;
    for (const std::size_t r : rows) sum += labels[r];
    out.base = sum / static_cast<double>(n);
    std::vector<double> residual(n);
    for (std::size_t i = 0; i < n; ++i) residual[i] = labels[rows[i]] - out.base;

    struct Pending {
        std::size_t node, begin, end, depth;
    };
    std::vector<std::size_t> index(n);
    std::vector<double> hs;
    std::vector<std::size_t> hc;
    for (std::size_t t = 0; t < o.trees; ++t) {
        std::vector<results::TreeNode> tree(1);
        std::iota(index.begin(), index.end(), std::size_t{0});
        // Each node owns a contiguous range of `index`; grown depth-first, the left child first.
        std::vector<Pending> stack{{0, 0, n, 0}};
        while (!stack.empty()) {
            const Pending p = stack.back();
            stack.pop_back();
            const std::size_t count = p.end - p.begin;
            double s = 0.0;
            for (std::size_t i = p.begin; i < p.end; ++i) s += residual[index[i]];
            std::size_t best_k = F, best_b = 0;
            double best_gain = 0.0;
            if (p.depth < o.depth && count >= 2 * o.min_leaf) {
                const double parent = s * s / static_cast<double>(count);
                for (std::size_t k = 0; k < F; ++k) {
                    const std::size_t bins = cuts[k].size() + 1;
                    if (bins < 2) continue;
                    hs.assign(bins, 0.0);
                    hc.assign(bins, 0);
                    for (std::size_t i = p.begin; i < p.end; ++i) {
                        const std::size_t r = index[i];
                        hs[bin[k][r]] += residual[r];
                        ++hc[bin[k][r]];
                    }
                    double sl = 0.0;
                    std::size_t cl = 0;
                    for (std::size_t b = 0; b + 1 < bins; ++b) {
                        sl += hs[b];
                        cl += hc[b];
                        const std::size_t cr = count - cl;
                        if (cl < o.min_leaf || cr < o.min_leaf) continue;
                        const double sr = s - sl;
                        const double gain =
                            sl * sl / static_cast<double>(cl) + sr * sr / static_cast<double>(cr) - parent;
                        if (gain > best_gain) {
                            best_gain = gain;
                            best_k = k;
                            best_b = b;
                        }
                    }
                }
            }
            if (best_k == F) {
                tree[p.node].value = o.shrinkage * s / static_cast<double>(count);
                continue;
            }
            const std::size_t k = best_k;
            out.importance[k] += best_gain;
            const auto first = index.begin() + static_cast<std::ptrdiff_t>(p.begin);
            const auto last = index.begin() + static_cast<std::ptrdiff_t>(p.end);
            const auto mid = std::stable_partition(first, last, [&](std::size_t r) { return bin[k][r] <= best_b; });
            const auto split = static_cast<std::size_t>(mid - index.begin());
            const std::size_t left = tree.size(), right = left + 1;
            tree.resize(tree.size() + 2);
            tree[p.node].feature = static_cast<std::int64_t>(k);
            tree[p.node].threshold = cuts[k][best_b];
            tree[p.node].left = static_cast<std::int64_t>(left);
            tree[p.node].right = static_cast<std::int64_t>(right);
            stack.push_back({right, split, p.end, p.depth + 1});
            stack.push_back({left, p.begin, split, p.depth + 1});
        }
        for (std::size_t i = 0; i < n; ++i) residual[i] -= tree[leaf_of(tree, &features[rows[i] * F])].value;
        out.trees.push_back(std::move(tree));
    }
    return out;
}

double predict_trees(double base, std::span<const std::vector<results::TreeNode>> trees, const double* row) {
    double y = base;
    for (const auto& tree : trees) y += tree[leaf_of(tree, row)].value;
    return y;
}

std::expected<Ridge, base::Error> fit_ridge(std::span<const double> features, std::span<const double> labels,
                                            std::span<const std::size_t> rows, double lambda) {
    if (rows.empty()) return std::unexpected(invalid("ridge regression needs samples"));
    Ridge r;
    r.mean.assign(F, 0.0);
    r.scale.assign(F, 0.0);
    r.categories.assign(F, {});
    const auto n = static_cast<double>(rows.size());
    for (std::size_t k = 0; k < F; ++k) {
        if (results::learn_feature_is_category(k)) {
            for (const std::size_t i : rows) r.categories[k].push_back(features[i * F + k]);
            std::ranges::sort(r.categories[k]);
            r.categories[k].erase(std::ranges::unique(r.categories[k]).begin(), r.categories[k].end());
            continue;
        }
        double s = 0.0;
        for (const std::size_t i : rows) s += features[i * F + k];
        r.mean[k] = s / n;
        double v = 0.0;
        for (const std::size_t i : rows) {
            const double d = features[i * F + k] - r.mean[k];
            v += d * d;
        }
        r.scale[k] = std::sqrt(v / n);
    }
    std::vector<double> x;
    ridge_row(r, &features[rows.front() * F], x);
    const std::size_t m = x.size();
    std::vector<double> A(m * m, 0.0), b(m, 0.0);
    for (const std::size_t i : rows) {
        ridge_row(r, &features[i * F], x);
        for (std::size_t p = 0; p < m; ++p) {
            b[p] += x[p] * labels[i];
            for (std::size_t q = 0; q <= p; ++q) A[p * m + q] += x[p] * x[q];
        }
    }
    for (std::size_t p = 1; p < m; ++p) A[p * m + p] += lambda;
    // Cholesky, the lower triangle in place.
    for (std::size_t j = 0; j < m; ++j) {
        double d = A[j * m + j];
        for (std::size_t k = 0; k < j; ++k) d -= A[j * m + k] * A[j * m + k];
        if (!(d > 0.0)) return std::unexpected(invalid("the ridge system is not positive definite"));
        A[j * m + j] = std::sqrt(d);
        for (std::size_t i = j + 1; i < m; ++i) {
            double s = A[i * m + j];
            for (std::size_t k = 0; k < j; ++k) s -= A[i * m + k] * A[j * m + k];
            A[i * m + j] = s / A[j * m + j];
        }
    }
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t k = 0; k < i; ++k) b[i] -= A[i * m + k] * b[k];
        b[i] /= A[i * m + i];
    }
    for (std::size_t i = m; i-- > 0;) {
        for (std::size_t k = i + 1; k < m; ++k) b[i] -= A[k * m + i] * b[k];
        b[i] /= A[i * m + i];
    }
    r.coefficients = std::move(b);
    return r;
}

double predict_ridge(const Ridge& r, const double* row) {
    std::vector<double> x;
    ridge_row(r, row, x);
    double y = 0.0;
    for (std::size_t p = 0; p < x.size(); ++p) y += r.coefficients[p] * x[p];
    return y;
}

double rank_correlation(std::span<const double> a, std::span<const double> b) {
    NITCAD_EXPECTS(a.size() == b.size());
    if (a.size() < 2) return 0.0;
    const std::vector<double> ra = ranks_of(a), rb = ranks_of(b);
    const auto n = static_cast<double>(a.size());
    const double mean = (n + 1.0) / 2.0;
    double sab = 0.0, saa = 0.0, sbb = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        sab += (ra[i] - mean) * (rb[i] - mean);
        saa += (ra[i] - mean) * (ra[i] - mean);
        sbb += (rb[i] - mean) * (rb[i] - mean);
    }
    if (!(saa > 0.0 && sbb > 0.0)) return 0.0;
    return sab / std::sqrt(saa * sbb);
}

std::vector<std::size_t> default_held_out(std::size_t devices) {
    std::vector<std::size_t> out;
    if (devices < 2) return out;
    for (std::size_t i = 3; i < devices; i += 4) out.push_back(i);
    if (out.empty()) out.push_back(devices - 1);
    return out;
}

std::expected<results::LearnedModel, base::Error> train_predictor(const results::TrainingSet& set,
                                                                  std::span<const std::size_t> held_out,
                                                                  const TreeOptions& o) {
    if (set.schema != results::learn_schema_version) {
        return std::unexpected(invalid("the training set's schema " + std::to_string(set.schema) +
                                       " is not this build's " + std::to_string(results::learn_schema_version)));
    }
    const std::size_t nd = set.devices.size();
    std::vector<char> held(nd, 0);
    for (const std::size_t d : held_out) {
        if (d >= nd) return std::unexpected(invalid("a held-out device is not in the manifest"));
        if (held[d] != 0) return std::unexpected(invalid("a held-out device is named twice"));
        held[d] = 1;
    }
    if (held_out.empty()) return std::unexpected(invalid("training needs at least one held-out device"));
    if (held_out.size() == nd) return std::unexpected(invalid("training needs at least one training device"));
    std::vector<std::size_t> train_rows, test_rows;
    for (std::size_t i = 0; i < set.rows; ++i) (held[set.device[i]] != 0 ? test_rows : train_rows).push_back(i);
    if (train_rows.empty()) return std::unexpected(invalid("the training devices have no samples"));
    if (test_rows.empty()) return std::unexpected(invalid("the held-out devices have no samples"));

    results::LearnedModel m;
    m.tree_count = o.trees;
    m.depth = o.depth;
    m.min_leaf = o.min_leaf;
    m.cut_points = o.cut_points;
    m.shrinkage = o.shrinkage;
    m.materials = set.materials;
    BoostedTrees bt = fit_boosted_trees(set.features, set.labels, train_rows, o);
    m.base = bt.base;
    m.trees = std::move(bt.trees);
    m.importance = std::move(bt.importance);
    m.feature_min.assign(F, std::numeric_limits<double>::infinity());
    m.feature_max.assign(F, -std::numeric_limits<double>::infinity());
    m.categories.assign(F, {});
    for (const std::size_t i : train_rows) {
        for (std::size_t k = 0; k < F; ++k) {
            const double v = set.features[i * F + k];
            m.feature_min[k] = std::min(m.feature_min[k], v);
            m.feature_max[k] = std::max(m.feature_max[k], v);
            if (results::learn_feature_is_category(k)) m.categories[k].push_back(v);
        }
    }
    for (std::size_t k = 0; k < F; ++k) {
        std::ranges::sort(m.categories[k]);
        m.categories[k].erase(std::ranges::unique(m.categories[k]).begin(), m.categories[k].end());
    }
    m.manifest_digest = manifest_digest_of(set.devices);
    for (std::size_t d = 0; d < nd; ++d) (held[d] != 0 ? m.held_out : m.training).push_back(set.devices[d]);
    m.training_samples = train_rows.size();
    // The build's facts only (LB13): the run's (start time, machine) would make two trainings on the
    // same data differ.
    {
        const results::Provenance p = provenance();
        m.provenance.version = p.version;
        m.provenance.revision = p.revision;
        m.provenance.dirty = p.dirty;
        m.provenance.build_type = p.build_type;
        m.provenance.compiler = p.compiler;
        m.provenance.eigen_version = p.eigen_version;
    }

    // The held-out metrics (LB10).
    auto ridge = fit_ridge(set.features, set.labels, train_rows);
    if (!ridge) return std::unexpected(ridge.error());
    std::vector<double> truth, trees, linear;
    std::size_t tp = 0, fp = 0, fn = 0, first = 0;
    for (const std::size_t i : test_rows) {
        const double* row = &set.features[i * F];
        const double y = predict_trees(m.base, m.trees, row);
        truth.push_back(set.labels[i]);
        trees.push_back(y);
        linear.push_back(predict_ridge(*ridge, row));
        if (set.pass[i] != 0) continue;
        ++first;
        const bool want = set.labels[i] >= 0.5, got = y >= 0.5;
        tp += want && got ? 1 : 0;
        fp += !want && got ? 1 : 0;
        fn += want && !got ? 1 : 0;
    }
    m.metrics.samples = test_rows.size();
    m.metrics.first_pass_samples = first;
    m.metrics.rank_correlation = rank_correlation(trees, truth);
    m.metrics.ridge_rank_correlation = rank_correlation(linear, truth);
    m.metrics.recall = tp + fn > 0 ? static_cast<double>(tp) / static_cast<double>(tp + fn) : 0.0;
    m.metrics.precision = tp + fp > 0 ? static_cast<double>(tp) / static_cast<double>(tp + fp) : 0.0;
    return m;
}

double predict_levels(const results::LearnedModel& model, const double* row) {
    return predict_trees(model.base, model.trees, row);
}

bool out_of_range(const results::LearnedModel& m, const double* row) {
    for (std::size_t k = 0; k < F; ++k) {
        const double v = row[k];
        if (results::learn_feature_is_category(k)) {
            if (!std::ranges::binary_search(m.categories[k], v)) return true;
            continue;
        }
        const double margin = 0.05 * (m.feature_max[k] - m.feature_min[k]);
        if (!(v >= m.feature_min[k] - margin && v <= m.feature_max[k] + margin)) return true;
    }
    return false;
}

double out_of_range_fraction(const results::LearnedModel& m, const EdgeFeatures& edges) {
    if (edges.rows == 0) return 0.0;
    std::size_t out = 0;
    for (std::size_t e = 0; e < edges.rows; ++e) out += out_of_range(m, &edges.features[e * F]) ? 1 : 0;
    return static_cast<double>(out) / static_cast<double>(edges.rows);
}

std::expected<SizeTargets, base::Error> learned_targets(const results::LearnedModel& m, const device::Device& d,
                                                        std::span<const std::vector<double>> points,
                                                        const BiasOptions& options) {
    auto f = edge_features(d, points, options, m.materials);
    if (!f) return std::unexpected(f.error());
    const std::size_t nn = d.mesh().node_count();
    SizeTargets t;
    t.size.assign(nn, std::numeric_limits<double>::infinity());
    t.normal.assign(nn, std::array<double, 3>{});
    const auto& edges = d.mesh().edges();
    for (std::size_t e = 0; e < f->rows; ++e) {
        const double l = predict_levels(m, &f->features[e * F]);
        if (!(l >= 0.5)) continue;
        const double target = f->length[e] * std::exp2(-l);
        for (const mesh::NodeId v : {edges[e].first, edges[e].second}) {
            const auto i = static_cast<std::size_t>(v);
            if (target < t.size[i]) {
                t.size[i] = target;
                t.normal[i] = f->direction[e];
            }
        }
    }
    return t;
}

}  // namespace NiTCAD::solve
