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
                       base::ThreadPool* pool)
    : pool_(pool),
      amg_(std::make_unique<AmgBackend>(amg, pool)),
      ilut_(std::make_unique<IlutBackend<double>>(drop_tolerance, fill_factor)) {}

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
    if (auto analyzed = ilut_->analyze(a); !analyzed) return analyzed;

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
    auto ilut = ilut_->factorize(values_);
    if (!ilut) return std::unexpected(std::move(ilut.error()));
    ++counts_.factorizations;
    return Factorization{.pivots = std::nullopt,
                         .perturbed_pivots = std::nullopt,
                         .fill = reduced_.nonzeros() + hierarchy->fill.value_or(0) +
                                 ilut->fill.value_or(0)};
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
    ilut_->solve(r2_, z2_);
    base::for_chunks(pool_, n, amg::vector_chunk, [&](std::size_t begin, std::size_t end) {
        for (std::size_t u = begin; u < end; ++u) x[u] = z1_[u] + z2_[u];
    });
    ++counts_.solves;
}

}  // namespace NiTCAD::linalg::detail
