#include "cpr.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>

#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/base/thread_pool.hpp"
#include "eigen_ilut.hpp"

namespace NiTCAD::linalg::detail {

namespace {

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = index, .value = std::nullopt}};
}

// Solves the m x m system a x = rhs (row-major a, overwritten; x returned in rhs) by Gaussian
// elimination with partial pivoting, the first largest pivot on a tie. False on an exactly zero
// pivot or a non-finite result.
bool solve_dense(std::span<double> a, std::span<double> rhs, std::size_t m) {
    for (std::size_t k = 0; k < m; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < m; ++i) {
            if (std::abs(a[i * m + k]) > std::abs(a[pivot * m + k])) pivot = i;
        }
        if (a[pivot * m + k] == 0.0) return false;
        if (pivot != k) {
            for (std::size_t j = 0; j < m; ++j) std::swap(a[k * m + j], a[pivot * m + j]);
            std::swap(rhs[k], rhs[pivot]);
        }
        for (std::size_t i = k + 1; i < m; ++i) {
            const double f = a[i * m + k] / a[k * m + k];
            for (std::size_t j = k; j < m; ++j) a[i * m + j] -= f * a[k * m + j];
            rhs[i] -= f * rhs[k];
        }
    }
    for (std::size_t k = m; k-- > 0;) {
        double s = rhs[k];
        for (std::size_t j = k + 1; j < m; ++j) s -= a[k * m + j] * rhs[j];
        rhs[k] = s / a[k * m + k];
        if (!std::isfinite(rhs[k])) return false;
    }
    return true;
}

}  // namespace

CprBackend::CprBackend(const AmgSettings& amg, double drop_tolerance, int fill_factor,
                       base::ThreadPool* pool, Stage2 stage2)
    : pool_(pool), stage2_(stage2), block_settings_(amg), amg_(std::make_unique<AmgBackend>(amg, pool)) {
    // block_amg keeps ILUT for retry_with_ilut.
    ilut_ = std::make_unique<IlutBackend<double>>(drop_tolerance, fill_factor);
    block_settings_.strength_decay = block_strength_decay;
}

CprBackend::~CprBackend() = default;

bool CprBackend::set_fields(std::span<const std::int32_t> labels) {
    if (std::ranges::equal(labels, labels_)) return false;
    labels_.assign(labels.begin(), labels.end());
    analyzed_ = false;
    return true;
}

std::expected<void, base::Error> CprBackend::analyze(const SparseMatrix& a) {
    NITCAD_EXPECTS(a.rows() == a.cols() && a.rows() > 0);
    analyzed_ = false;
    const auto n = static_cast<std::size_t>(a.rows());
    if (labels_.size() != n) {
        return std::unexpected(invalid("the cpr preconditioner needs one field label per unknown "
                                       "(set_field_labels)",
                                       labels_.size()));
    }
    // Nodes from the labels.
    potential_.clear();
    node_end_.clear();
    node_of_.assign(n, -1);
    for (std::size_t u = 0; u < n; ++u) {
        if (labels_[u] == 0) {
            node_of_[u] = static_cast<Index>(potential_.size());
            potential_.push_back(u);
            node_end_.push_back(u + 1);
        } else if (labels_[u] > 0) {
            if (potential_.empty() || node_end_.back() != u) {
                return std::unexpected(
                    invalid("cpr: a coupled field's unknown has no potential before it in its node", u));
            }
            node_of_[u] = static_cast<Index>(potential_.size() - 1);
            node_end_.back() = u + 1;
        }
    }
    if (potential_.empty()) {
        return std::unexpected(invalid("cpr needs at least one potential unknown"));
    }
    if (ilut_) {
        if (auto analyzed = ilut_->analyze(a); !analyzed) return analyzed;
    }

    offsets_.assign(a.row_offsets().begin(), a.row_offsets().end());
    columns_.assign(a.col_indices().begin(), a.col_indices().end());
    values_.assign(a.nonzeros(), 0.0);

    // A_p's pattern: potential row g holds the nodes of its row's columns.
    const std::size_t np = potential_.size();
    std::vector<Triplet> entries;
    for (std::size_t g = 0; g < np; ++g) {
        const std::size_t r = potential_[g];
        for (auto k = static_cast<std::size_t>(offsets_[r]); k < static_cast<std::size_t>(offsets_[r + 1]);
             ++k) {
            const Index node = node_of_[static_cast<std::size_t>(columns_[k])];
            if (node >= 0) entries.push_back({static_cast<Index>(g), node, 0.0});
        }
        entries.push_back({static_cast<Index>(g), static_cast<Index>(g), 0.0});  // the diagonal
    }
    auto reduced = SparseMatrix::from_triplets(static_cast<Index>(np), static_cast<Index>(np), entries);
    if (!reduced) return std::unexpected(std::move(reduced.error()));
    reduced_ = std::move(*reduced);
    target_.assign(a.nonzeros(), -1);
    const std::span<const Index> ro = reduced_.row_offsets();
    const std::span<const Index> rc = reduced_.col_indices();
    for (std::size_t g = 0; g < np; ++g) {
        const std::size_t r = potential_[g];
        for (auto k = static_cast<std::size_t>(offsets_[r]); k < static_cast<std::size_t>(offsets_[r + 1]);
             ++k) {
            const Index node = node_of_[static_cast<std::size_t>(columns_[k])];
            if (node < 0) continue;
            const auto first = rc.begin() + ro[g];
            const auto last = rc.begin() + ro[g + 1];
            const auto it = std::lower_bound(first, last, node);
            NITCAD_EXPECTS(it != last && *it == node);
            target_[k] = static_cast<Index>(it - rc.begin());
        }
    }
    if (auto analyzed = amg_->analyze(reduced_); !analyzed) return analyzed;

    response_.assign(n, 0.0);
    rows_p_.assign(np, 1.0);
    cols_p_.assign(np, 1.0);
    v_p_.assign(np, 0.0);
    y_p_.assign(np, 0.0);
    for (auto* w : {&z1_, &r2_, &z2_}) w->assign(n, 0.0);

    // block_amg: one block per field label, its pattern, and every entry's place in it.
    block_labels_.clear();
    blocks_.clear();
    outside_.clear();
    block_entry_.assign(a.nonzeros(), -1);
    block_of_.assign(n, -1);
    position_.assign(n, 0);
    if (stage2_ == Stage2::block_amg) {
        for (const std::int32_t l : labels_) {
            if (l >= 0) block_labels_.push_back(l);
        }
        std::ranges::sort(block_labels_);
        const auto [last, end] = std::ranges::unique(block_labels_);
        block_labels_.erase(last, end);
        blocks_.resize(block_labels_.size());
        for (std::size_t u = 0; u < n; ++u) {
            if (labels_[u] < 0) {
                outside_.push_back(u);
                continue;
            }
            const auto b = static_cast<std::size_t>(std::ranges::lower_bound(block_labels_, labels_[u]) -
                                                    block_labels_.begin());
            block_of_[u] = static_cast<Index>(b);
            position_[u] = blocks_[b].unknowns.size();
            blocks_[b].unknowns.push_back(u);
        }
        for (std::size_t b = 0; b < blocks_.size(); ++b) {
            Block& block = blocks_[b];
            std::vector<Triplet> block_entries;
            for (const std::size_t r : block.unknowns) {
                for (auto k = static_cast<std::size_t>(offsets_[r]); k < static_cast<std::size_t>(offsets_[r + 1]);
                     ++k) {
                    const auto c = static_cast<std::size_t>(columns_[k]);
                    if (block_of_[c] == static_cast<Index>(b)) {
                        block_entries.push_back({static_cast<Index>(position_[r]), static_cast<Index>(position_[c]), 0.0});
                    }
                }
                block_entries.push_back({static_cast<Index>(position_[r]), static_cast<Index>(position_[r]), 0.0});
            }
            const auto m = static_cast<Index>(block.unknowns.size());
            auto matrix = SparseMatrix::from_triplets(m, m, block_entries);
            if (!matrix) return std::unexpected(std::move(matrix.error()));
            block.matrix = std::move(*matrix);
            const std::span<const Index> bo = block.matrix.row_offsets();
            const std::span<const Index> bc = block.matrix.col_indices();
            for (const std::size_t r : block.unknowns) {
                const std::size_t pr = position_[r];
                for (auto k = static_cast<std::size_t>(offsets_[r]); k < static_cast<std::size_t>(offsets_[r + 1]);
                     ++k) {
                    const auto c = static_cast<std::size_t>(columns_[k]);
                    if (block_of_[c] != static_cast<Index>(b)) continue;
                    const auto first = bc.begin() + bo[pr];
                    const auto it = std::lower_bound(first, bc.begin() + bo[pr + 1], static_cast<Index>(position_[c]));
                    block_entry_[k] = static_cast<Index>(it - bc.begin());
                }
            }
            block.amg = std::make_unique<SymmetrizedAmg>(block_settings_, pool_);
            if (auto analyzed = block.amg->analyze(block.matrix); !analyzed) return analyzed;
            block.rhs.assign(block.unknowns.size(), 0.0);
            block.z.assign(block.unknowns.size(), 0.0);
        }
    }
    analyzed_ = true;
    ++counts_.analyses;
    return {};
}

void CprBackend::set_scales(std::span<const double> rows, std::span<const double> columns,
                            bool applied) {
    NITCAD_EXPECTS(analyzed_ && rows.size() == response_.size() && columns.size() == rows.size());
    for (std::size_t g = 0; g < potential_.size(); ++g) {
        rows_p_[g] = rows[potential_[g]];
        cols_p_[g] = columns[potential_[g]];
    }
    scaled_ = applied;
}

std::expected<Factorization, base::Error> CprBackend::factorize(std::span<const double> csr_values) {
    NITCAD_EXPECTS(analyzed_ && csr_values.size() == values_.size());
    std::ranges::copy(csr_values, values_.begin());

    // The responses e of every node's other unknowns (serial, node order).
    const auto entry = [&](std::size_t row, std::size_t col) {
        for (auto k = static_cast<std::size_t>(offsets_[row]); k < static_cast<std::size_t>(offsets_[row + 1]);
             ++k) {
            if (static_cast<std::size_t>(columns_[k]) == col) return values_[k];
        }
        return 0.0;
    };
    for (std::size_t g = 0; g < potential_.size(); ++g) {
        const std::size_t s = potential_[g];
        response_[s] = 1.0;
        const std::size_t m = node_end_[g] - s - 1;
        if (m == 0) continue;
        block_.assign(m * m + m, 0.0);
        const std::span<double> b(block_.data(), m * m);
        const std::span<double> e(block_.data() + m * m, m);
        for (std::size_t i = 0; i < m; ++i) {
            for (std::size_t j = 0; j < m; ++j) b[i * m + j] = entry(s + 1 + i, s + 1 + j);
            e[i] = -entry(s + 1 + i, s);
        }
        const bool solved = solve_dense(b, e, m);
        for (std::size_t i = 0; i < m; ++i) response_[s + 1 + i] = solved ? e[i] : 0.0;
    }

    // A_p = R A P.
    const std::span<double> reduced = reduced_.values();
    std::ranges::fill(reduced, 0.0);
    for (std::size_t g = 0; g < potential_.size(); ++g) {
        const std::size_t r = potential_[g];
        for (auto k = static_cast<std::size_t>(offsets_[r]); k < static_cast<std::size_t>(offsets_[r + 1]);
             ++k) {
            if (target_[k] < 0) continue;
            reduced[static_cast<std::size_t>(target_[k])] +=
                values_[k] * response_[static_cast<std::size_t>(columns_[k])];
        }
    }
    for (const double v : reduced) {
        if (!std::isfinite(v)) {
            return std::unexpected(base::Error{base::ErrorCode::inaccurate_solve,
                                               "cpr: the potential system is not finite",
                                               std::nullopt});
        }
    }

    amg_->set_scales(rows_p_, cols_p_, scaled_);
    auto hierarchy = amg_->factorize(reduced);
    if (!hierarchy) return std::unexpected(std::move(hierarchy.error()));
    std::size_t stage2_fill = 0;
    // block_amg: the swept blocks' AMG; ILUT only on retry_with_ilut.
    use_ilut_ = stage2_ == Stage2::ilut;
    for (std::size_t b = 0; b < blocks_.size() && !use_ilut_; ++b) {
        if (sweep_ != Sweep::all && block_labels_[b] == 0) continue;  // not swept
        Block& block = blocks_[b];
        const std::span<double> v = block.matrix.values();
        std::ranges::fill(v, 0.0);
        for (const std::size_t r : block.unknowns) {
            for (auto k = static_cast<std::size_t>(offsets_[r]); k < static_cast<std::size_t>(offsets_[r + 1]); ++k) {
                if (block_entry_[k] >= 0) v[static_cast<std::size_t>(block_entry_[k])] = values_[k];
            }
        }
        auto built = block.amg->factorize(v);
        if (!built) return std::unexpected(std::move(built.error()));
        stage2_fill += block.matrix.nonzeros() + built->fill.value_or(0);
    }
    if (use_ilut_) {
        auto ilut = ilut_->factorize(values_);
        if (!ilut) return std::unexpected(std::move(ilut.error()));
        stage2_fill = ilut->fill.value_or(0);
    }
    ++counts_.factorizations;
    return Factorization{.pivots = std::nullopt,
                         .perturbed_pivots = std::nullopt,
                         .fill = reduced_.nonzeros() + hierarchy->fill.value_or(0) +
                                 stage2_fill};
}

std::expected<bool, base::Error> CprBackend::retry_with_ilut() {
    if (use_ilut_) return false;
    auto ilut = ilut_->factorize(values_);
    if (!ilut) return std::unexpected(std::move(ilut.error()));
    use_ilut_ = true;
    ++counts_.stage2_retries;
    return true;
}

void CprBackend::solve(std::span<const double> b, std::span<double> x) {
    NITCAD_EXPECTS(analyzed_ && b.size() == response_.size() && x.size() == b.size());
    const std::size_t n = b.size();
    // Stage 1: AMG on the potential system, prolonged with the responses.
    for (std::size_t g = 0; g < potential_.size(); ++g) v_p_[g] = b[potential_[g]];
    amg_->solve(v_p_, y_p_);
    base::for_chunks(pool_, n, amg::vector_chunk, [&](std::size_t begin, std::size_t end) {
        for (std::size_t u = begin; u < end; ++u) {
            const Index node = node_of_[u];
            z1_[u] = node < 0 ? 0.0 : response_[u] * y_p_[static_cast<std::size_t>(node)];
        }
    });
    // Stage 2: ILUT on the residual stage 1 leaves.
    base::for_chunks(pool_, n, amg::row_chunk, [&](std::size_t begin, std::size_t end) {
        for (std::size_t r = begin; r < end; ++r) {
            double s = 0.0;
            for (auto k = static_cast<std::size_t>(offsets_[r]); k < static_cast<std::size_t>(offsets_[r + 1]);
                 ++k) {
                s += values_[k] * z1_[static_cast<std::size_t>(columns_[k])];
            }
            r2_[r] = b[r] - s;
        }
    });
    if (use_ilut_) {
        ilut_->solve(r2_, z2_);
    } else {
        solve_blocks();
    }
    base::for_chunks(pool_, n, amg::vector_chunk, [&](std::size_t begin, std::size_t end) {
        for (std::size_t u = begin; u < end; ++u) x[u] = z1_[u] + z2_[u];
    });
    ++counts_.solves;
}

void CprBackend::solve_blocks() {
    std::ranges::fill(z2_, 0.0);
    for (std::size_t b = 0; b < blocks_.size(); ++b) {
        Block& block = blocks_[b];
        if (sweep_ != Sweep::all && block_labels_[b] == 0) continue;  // stage 1 owns the potential
        // The field's right-hand side: r2 less the couplings to the fields already swept (every
        // other field's z2 is still zero).
        base::for_chunks(pool_, block.unknowns.size(), amg::row_chunk, [&](std::size_t begin, std::size_t end) {
            for (std::size_t i = begin; i < end; ++i) {
                const std::size_t r = block.unknowns[i];
                double s = r2_[r];
                for (auto k = static_cast<std::size_t>(offsets_[r]); k < static_cast<std::size_t>(offsets_[r + 1]);
                     ++k) {
                    if (block_entry_[k] < 0) s -= values_[k] * z2_[static_cast<std::size_t>(columns_[k])];
                }
                block.rhs[i] = s;
            }
        });
        block.amg->solve(block.rhs, block.z);
        base::for_chunks(pool_, block.unknowns.size(), amg::vector_chunk, [&](std::size_t begin, std::size_t end) {
            for (std::size_t i = begin; i < end; ++i) z2_[block.unknowns[i]] = block.z[i];
        });
    }
    // The unknowns in no node: one Gauss-Seidel step each, in order.
    for (const std::size_t u : outside_) {
        double s = r2_[u];
        double diagonal = 0.0;
        for (auto k = static_cast<std::size_t>(offsets_[u]); k < static_cast<std::size_t>(offsets_[u + 1]); ++k) {
            const auto c = static_cast<std::size_t>(columns_[k]);
            if (c == u) {
                diagonal = values_[k];
            } else {
                s -= values_[k] * z2_[c];
            }
        }
        z2_[u] = diagonal != 0.0 ? s / diagonal : 0.0;
    }
}

}  // namespace NiTCAD::linalg::detail
