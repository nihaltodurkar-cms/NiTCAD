#include "amg.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>

#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/base/thread_pool.hpp"
#include "eigen_sparse_lu.hpp"

namespace NiTCAD::linalg::detail {

namespace {

base::Error build_error(const char* message) {
    return base::Error{base::ErrorCode::inaccurate_solve, message, std::nullopt};
}

std::size_t at(Index i) noexcept { return static_cast<std::size_t>(i); }

// 1 / a_ii; a row without a diagonal value takes its l1 norm (sign +), an all-zero row 0.
std::vector<double> inverse_diagonal(const Csr& a) {
    std::vector<double> inv(a.rows, 0.0);
    for (std::size_t i = 0; i < a.rows; ++i) {
        double sum = 0.0;
        double diagonal = 0.0;
        for (std::size_t k = at(a.offsets[i]); k < at(a.offsets[i + 1]); ++k) {
            sum += std::abs(a.values[k]);
            if (at(a.columns[k]) == i) diagonal = a.values[k];
        }
        if (diagonal != 0.0) {
            inv[i] = 1.0 / diagonal;
        } else if (sum > 0.0) {
            inv[i] = 1.0 / sum;
        }
    }
    return inv;
}

// The largest eigenvalue modulus of D^-1 A, by power iteration from a fixed vector (serial, so the
// same bits every time), times a safety factor; 1 if the estimate is not a positive finite number.
double spectral_radius(const Csr& a, std::span<const double> inv_d) {
    const std::size_t n = a.rows;
    std::vector<double> v(n);
    std::vector<double> w(n);
    for (std::size_t i = 0; i < n; ++i) {
        // A fixed pseudo-random start in [-1, 1) (splitmix64 of the index): every mode present.
        std::uint64_t z = (static_cast<std::uint64_t>(i) + 1) * 0x9E3779B97F4A7C15ull;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        z ^= z >> 31;
        v[i] = static_cast<double>(z >> 11) * 0x1.0p-52 - 1.0;
    }
    const auto norm = [](std::span<const double> u) {
        double s = 0.0;
        for (const double x : u) s += x * x;
        return std::sqrt(s);
    };
    double lambda = 0.0;
    double v_norm = norm(v);
    for (int step = 0; step < amg::power_iterations && v_norm > 0.0; ++step) {
        amg::multiply(nullptr, a, v, w);
        for (std::size_t i = 0; i < n; ++i) w[i] *= inv_d[i];
        const double w_norm = norm(w);
        lambda = w_norm / v_norm;
        for (std::size_t i = 0; i < n; ++i) v[i] = w_norm > 0.0 ? w[i] / w_norm : 0.0;
        v_norm = w_norm > 0.0 ? 1.0 : 0.0;
    }
    lambda *= amg::spectral_safety;
    return lambda > 0.0 && std::isfinite(lambda) ? lambda : 1.0;
}

bool all_finite(std::span<const double> v) {
    return std::ranges::all_of(v, [](double x) { return std::isfinite(x); });
}

}  // namespace

namespace amg {

void multiply(base::ThreadPool* pool, const Csr& a, std::span<const double> x, std::span<double> y) {
    NITCAD_EXPECTS(x.size() == a.cols && y.size() == a.rows);
    base::for_chunks(pool, a.rows, row_chunk, [&](std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            double s = 0.0;
            for (std::size_t k = at(a.offsets[i]); k < at(a.offsets[i + 1]); ++k) {
                s += a.values[k] * x[at(a.columns[k])];
            }
            y[i] = s;
        }
    });
}

Csr transpose(const Csr& a) {
    Csr t;
    t.rows = a.cols;
    t.cols = a.rows;
    t.offsets.assign(a.cols + 1, 0);
    for (const Index c : a.columns) ++t.offsets[at(c) + 1];
    for (std::size_t c = 0; c < a.cols; ++c) t.offsets[c + 1] += t.offsets[c];
    t.columns.resize(a.nonzeros());
    t.values.resize(a.nonzeros());
    std::vector<Index> next(t.offsets.begin(), t.offsets.end() - 1);
    // Rows of A in order, so each row of A^T receives its columns in increasing order.
    for (std::size_t i = 0; i < a.rows; ++i) {
        for (std::size_t k = at(a.offsets[i]); k < at(a.offsets[i + 1]); ++k) {
            const std::size_t pos = at(next[at(a.columns[k])]++);
            t.columns[pos] = static_cast<Index>(i);
            t.values[pos] = a.values[k];
        }
    }
    return t;
}

Csr product(const Csr& a, const Csr& b) {
    NITCAD_EXPECTS(a.cols == b.rows);
    Csr c;
    c.rows = a.rows;
    c.cols = b.cols;
    c.offsets.assign(a.rows + 1, 0);
    std::vector<std::size_t> marker(b.cols, std::numeric_limits<std::size_t>::max());
    std::vector<double> accumulator(b.cols, 0.0);
    std::vector<Index> row;
    for (std::size_t i = 0; i < a.rows; ++i) {
        row.clear();
        for (std::size_t k = at(a.offsets[i]); k < at(a.offsets[i + 1]); ++k) {
            const double aik = a.values[k];
            const std::size_t kk = at(a.columns[k]);
            for (std::size_t q = at(b.offsets[kk]); q < at(b.offsets[kk + 1]); ++q) {
                const std::size_t j = at(b.columns[q]);
                if (marker[j] != i) {
                    marker[j] = i;
                    accumulator[j] = aik * b.values[q];
                    row.push_back(static_cast<Index>(j));
                } else {
                    accumulator[j] += aik * b.values[q];
                }
            }
        }
        std::ranges::sort(row);
        for (const Index j : row) {
            c.columns.push_back(j);
            c.values.push_back(accumulator[at(j)]);
        }
        c.offsets[i + 1] = static_cast<Index>(c.columns.size());
    }
    return c;
}

Coarsening coarsen(const Csr& a, double theta) {
    return coarsen(a, StrengthMeasure::scaled, theta, {});
}

Coarsening coarsen(const Csr& a, StrengthMeasure measure, double theta,
                   std::span<const double> candidate) {
    const std::size_t n = a.rows;
    Coarsening out;
    out.strong.assign(a.nonzeros(), 0);
    std::vector<double> diagonal(n, 0.0);
    std::vector<std::uint8_t> isolated(n, 1);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = at(a.offsets[i]); k < at(a.offsets[i + 1]); ++k) {
            if (at(a.columns[k]) == i) {
                diagonal[i] = std::abs(a.values[k]);
            } else if (a.values[k] != 0.0) {
                isolated[i] = 0;
            }
        }
    }
    // Directed strength, then the symmetric graph S (with each edge's larger magnitude).
    Csr directed;
    directed.rows = directed.cols = n;
    directed.offsets.assign(n + 1, 0);

    // The evolution measure's state (Unit 63 MC1 (2)): the transposed pattern (rows reaching a
    // support node), omega, and the work vectors (zero outside the current support).
    Csr reach;
    double omega = 0.0;
    std::vector<double> z, z_next, measure_of;
    std::vector<std::size_t> support, next;
    std::vector<std::size_t> stamp;
    std::size_t stamp_id = 0;
    if (measure == StrengthMeasure::evolution) {
        reach = transpose(a);
        double gershgorin = 0.0;
        for (std::size_t r = 0; r < n; ++r) {
            const double ar = diagonal[r];
            if (ar == 0.0) continue;
            double sum = 0.0;
            for (std::size_t k = at(a.offsets[r]); k < at(a.offsets[r + 1]); ++k) sum += std::abs(a.values[k]);
            gershgorin = std::max(gershgorin, sum / ar);
        }
        omega = gershgorin > 0.0 ? 1.0 / gershgorin : 0.0;
        z.assign(n, 0.0);
        z_next.assign(n, 0.0);
        stamp.assign(n, 0);
    }
    const auto b = [&](std::size_t u) { return candidate.empty() ? 1.0 : candidate[u]; };
    const auto signed_diagonal = [&](std::size_t r) {
        for (std::size_t k = at(a.offsets[r]); k < at(a.offsets[r + 1]); ++k) {
            if (at(a.columns[k]) == r) return a.values[k];
        }
        return 0.0;
    };
    // z = (I - omega D^-1 A)^steps e_i over a growing support; m_ij for row i's entries.
    const auto evolve = [&](std::size_t i) {
        support.assign(1, i);
        z[i] = 1.0;
        for (int step = 0; step < evolution_steps; ++step) {
            ++stamp_id;
            next.clear();
            for (const std::size_t s : support) {
                if (stamp[s] != stamp_id) { stamp[s] = stamp_id; next.push_back(s); }
                for (std::size_t k = at(reach.offsets[s]); k < at(reach.offsets[s + 1]); ++k) {
                    const std::size_t r = at(reach.columns[k]);
                    if (stamp[r] != stamp_id) { stamp[r] = stamp_id; next.push_back(r); }
                }
            }
            std::ranges::sort(next);
            for (const std::size_t r : next) {
                double az = 0.0;
                for (std::size_t k = at(a.offsets[r]); k < at(a.offsets[r + 1]); ++k) {
                    az += a.values[k] * z[at(a.columns[k])];
                }
                const double ar = signed_diagonal(r);
                z_next[r] = ar != 0.0 ? z[r] - omega * az / ar : z[r];
            }
            for (const std::size_t r : next) z[r] = z_next[r];
            support.swap(next);
        }
        measure_of.assign(at(a.offsets[i + 1]) - at(a.offsets[i]), std::numeric_limits<double>::infinity());
        double smallest = std::numeric_limits<double>::infinity();
        for (std::size_t k = at(a.offsets[i]); k < at(a.offsets[i + 1]); ++k) {
            const std::size_t j = at(a.columns[k]);
            if (j == i || isolated[j] != 0 || a.values[k] == 0.0 || z[i] == 0.0 || b(j) == 0.0) continue;
            const double m = std::abs(1.0 - (z[j] * b(i)) / (z[i] * b(j)));
            measure_of[k - at(a.offsets[i])] = m;
            smallest = std::min(smallest, m);
        }
        for (const std::size_t r : support) { z[r] = 0.0; z_next[r] = 0.0; }
        return std::max(smallest, evolution_floor);
    };

    for (std::size_t i = 0; i < n; ++i) {
        if (isolated[i] == 0) {
            double row_largest = 0.0;
            double evolution_min = 0.0;
            if (measure == StrengthMeasure::row_max) {
                for (std::size_t k = at(a.offsets[i]); k < at(a.offsets[i + 1]); ++k) {
                    const std::size_t j = at(a.columns[k]);
                    if (j != i && isolated[j] == 0) row_largest = std::max(row_largest, std::abs(a.values[k]));
                }
            } else if (measure == StrengthMeasure::evolution) {
                evolution_min = evolve(i);
            }
            for (std::size_t k = at(a.offsets[i]); k < at(a.offsets[i + 1]); ++k) {
                const std::size_t j = at(a.columns[k]);
                const double v = std::abs(a.values[k]);
                bool strong = false;
                switch (measure) {
                    case StrengthMeasure::scaled:
                        strong = v >= theta * std::sqrt(diagonal[i] * diagonal[j]);
                        break;
                    case StrengthMeasure::row_max:
                        strong = v >= theta * row_largest;
                        break;
                    case StrengthMeasure::evolution:
                        strong = measure_of[k - at(a.offsets[i])] <= theta * evolution_min;
                        break;
                }
                if (j != i && isolated[j] == 0 && v > 0.0 && strong) {
                    directed.columns.push_back(static_cast<Index>(j));
                    // The scale-free strength, so pass 2 compares like with like.
                    const double d = std::sqrt(diagonal[i] * diagonal[j]);
                    directed.values.push_back(d > 0.0 ? v / d : std::numeric_limits<double>::infinity());
                }
            }
        }
        directed.offsets[i + 1] = static_cast<Index>(directed.columns.size());
    }
    const Csr reverse = transpose(directed);
    Csr s;
    s.rows = s.cols = n;
    s.offsets.assign(n + 1, 0);
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t p = at(directed.offsets[i]);
        std::size_t q = at(reverse.offsets[i]);
        const std::size_t pe = at(directed.offsets[i + 1]);
        const std::size_t qe = at(reverse.offsets[i + 1]);
        while (p < pe || q < qe) {
            const Index cp = p < pe ? directed.columns[p] : std::numeric_limits<Index>::max();
            const Index cq = q < qe ? reverse.columns[q] : std::numeric_limits<Index>::max();
            if (cp == cq) {
                s.columns.push_back(cp);
                s.values.push_back(std::max(directed.values[p++], reverse.values[q++]));
            } else if (cp < cq) {
                s.columns.push_back(cp);
                s.values.push_back(directed.values[p++]);
            } else {
                s.columns.push_back(cq);
                s.values.push_back(reverse.values[q++]);
            }
        }
        s.offsets[i + 1] = static_cast<Index>(s.columns.size());
    }
    // The strong mask over A's entries: those in S (both sorted, so a merge).
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t q = at(s.offsets[i]);
        for (std::size_t k = at(a.offsets[i]); k < at(a.offsets[i + 1]); ++k) {
            while (q < at(s.offsets[i + 1]) && s.columns[q] < a.columns[k]) ++q;
            if (q < at(s.offsets[i + 1]) && s.columns[q] == a.columns[k]) out.strong[k] = 1;
        }
    }

    constexpr Index free_node = -1;
    std::vector<Index>& agg = out.aggregate;
    agg.assign(n, free_node);
    Index count = 0;
    // Pass 1: a node whose strong neighbours are all free starts an aggregate with them.
    for (std::size_t i = 0; i < n; ++i) {
        if (isolated[i] != 0 || agg[i] != free_node) continue;
        const std::size_t begin = at(s.offsets[i]);
        const std::size_t end = at(s.offsets[i + 1]);
        if (begin == end) continue;  // no strong neighbour: pass 3 makes it a singleton
        bool all_free = true;
        for (std::size_t k = begin; k < end && all_free; ++k) all_free = agg[at(s.columns[k])] == free_node;
        if (!all_free) continue;
        agg[i] = count;
        for (std::size_t k = begin; k < end; ++k) agg[at(s.columns[k])] = count;
        ++count;
    }
    // Pass 2: join the strongest neighbour's aggregate of pass 1.
    const std::vector<Index> first = agg;
    for (std::size_t i = 0; i < n; ++i) {
        if (isolated[i] != 0 || agg[i] != free_node) continue;
        double best = -1.0;
        Index chosen = free_node;
        for (std::size_t k = at(s.offsets[i]); k < at(s.offsets[i + 1]); ++k) {
            const Index target = first[at(s.columns[k])];
            if (target != free_node && s.values[k] > best) {
                best = s.values[k];
                chosen = target;
            }
        }
        agg[i] = chosen;
    }
    // Pass 3: what is still free starts an aggregate with its free strong neighbours (a node
    // without any is a singleton: excluding such nodes instead measured worse on coarse levels,
    // 16 against 12 iterations on the 48^3 Poisson problem).
    for (std::size_t i = 0; i < n; ++i) {
        if (isolated[i] != 0 || agg[i] != free_node) continue;
        agg[i] = count;
        for (std::size_t k = at(s.offsets[i]); k < at(s.offsets[i + 1]); ++k) {
            if (agg[at(s.columns[k])] == free_node) agg[at(s.columns[k])] = count;
        }
        ++count;
    }
    out.aggregates = static_cast<std::size_t>(count);
    return out;
}

namespace {

constexpr std::size_t no_entry = std::numeric_limits<std::size_t>::max();

// The position of entry (row, col) in a's arrays, or no_entry.
std::size_t find_entry(const Csr& a, std::size_t row, Index col) {
    const auto first = a.columns.begin() + a.offsets[row];
    const auto last = a.columns.begin() + a.offsets[row + 1];
    const auto it = std::lower_bound(first, last, col);
    return it != last && *it == col ? static_cast<std::size_t>(it - a.columns.begin()) : no_entry;
}

// Entry k = (i, j) of a and its partner t = (j, i) form a usable pair.
bool usable(const Csr& a, std::size_t k, std::size_t t) {
    if (t == no_entry) return false;
    const double x = a.values[k];
    const double y = a.values[t];
    return x != 0.0 && y != 0.0 && (x > 0.0) == (y > 0.0) && std::isfinite(x) && std::isfinite(y);
}

}  // namespace

Symmetrization symmetrize(const Csr& a) {
    NITCAD_EXPECTS(a.rows == a.cols);
    const std::size_t n = a.rows;
    Symmetrization s;
    std::vector<double> ln(n, 0.0);
    std::vector<std::uint8_t> seen(n, 0);
    std::vector<std::size_t> queue;
    queue.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = at(a.offsets[i]); k < at(a.offsets[i + 1]); ++k) {
            const std::size_t j = at(a.columns[k]);
            if (j != i && !usable(a, k, find_entry(a, j, static_cast<Index>(i)))) ++s.unpaired;
        }
    }
    for (std::size_t root = 0; root < n; ++root) {
        if (seen[root]) continue;
        ++s.components;
        seen[root] = 1;
        queue.clear();
        queue.push_back(root);
        for (std::size_t q = 0; q < queue.size(); ++q) {
            const std::size_t i = queue[q];
            for (std::size_t k = at(a.offsets[i]); k < at(a.offsets[i + 1]); ++k) {
                const std::size_t j = at(a.columns[k]);
                if (j == i || seen[j]) continue;
                const std::size_t t = find_entry(a, j, static_cast<Index>(i));
                if (!usable(a, k, t)) continue;
                ln[j] = ln[i] + 0.5 * (std::log(std::abs(a.values[k])) - std::log(std::abs(a.values[t])));
                seen[j] = 1;
                queue.push_back(j);
            }
        }
        double mean = 0.0;
        for (const std::size_t u : queue) mean += ln[u];
        mean /= static_cast<double>(queue.size());
        for (const std::size_t u : queue) {
            ln[u] -= mean;
            if (std::abs(ln[u]) > symmetrization_clamp) {
                ln[u] = std::copysign(symmetrization_clamp, ln[u]);
                ++s.clamped;
            }
        }
    }
    s.d.resize(n);
    for (std::size_t i = 0; i < n; ++i) s.d[i] = std::exp(ln[i]);
    return s;
}

std::vector<double> symmetric_part(const Csr& a, std::span<const double> d) {
    NITCAD_EXPECTS(a.rows == a.cols && d.size() == a.rows);
    std::vector<double> v(a.nonzeros());
    for (std::size_t i = 0; i < a.rows; ++i) {
        for (std::size_t k = at(a.offsets[i]); k < at(a.offsets[i + 1]); ++k) {
            const std::size_t j = at(a.columns[k]);
            const double s_ij = d[i] * a.values[k] / d[j];
            const std::size_t t = find_entry(a, j, static_cast<Index>(i));
            v[k] = j != i && usable(a, k, t) ? 0.5 * (s_ij + d[j] * a.values[t] / d[i]) : s_ij;
        }
    }
    return v;
}

}  // namespace amg

AmgBackend::AmgBackend(const AmgSettings& settings, base::ThreadPool* pool)
    : settings_(settings), pool_(pool) {}

AmgBackend::~AmgBackend() = default;

std::expected<void, base::Error> AmgBackend::analyze(const SparseMatrix& a) {
    NITCAD_EXPECTS(a.rows() == a.cols() && a.rows() > 0);
    built_ = false;
    coarsened_ = false;
    coarse_smoothed_ = false;
    coarse_analyzed_ = false;
    coarsening_.clear();
    levels_.clear();
    Level& fine = levels_.emplace_back();
    fine.a.rows = fine.a.cols = static_cast<std::size_t>(a.rows());
    fine.a.offsets.assign(a.row_offsets().begin(), a.row_offsets().end());
    fine.a.columns.assign(a.col_indices().begin(), a.col_indices().end());
    fine.a.values.assign(a.nonzeros(), 0.0);
    candidate_.assign(fine.a.rows, 1.0);
    similarity_.assign(fine.a.rows, 1.0);
    ++counts_.analyses;
    return {};
}

void AmgBackend::set_scales(std::span<const double> rows, std::span<const double> columns,
                            bool applied) {
    NITCAD_EXPECTS(!levels_.empty() && rows.size() == levels_.front().a.rows &&
                   columns.size() == rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) {
        similarity_[i] = applied ? std::sqrt(columns[i] / rows[i]) : 1.0;
        candidate_[i] = applied ? 1.0 / std::sqrt(rows[i] * columns[i]) : 1.0;
    }
}

std::vector<std::size_t> AmgBackend::level_sizes() const {
    std::vector<std::size_t> sizes;
    sizes.reserve(levels_.size());
    for (const Level& level : levels_) sizes.push_back(level.a.rows);
    return sizes;
}

std::expected<Factorization, base::Error> AmgBackend::factorize(std::span<const double> csr_values) {
    NITCAD_EXPECTS(!levels_.empty());
    NITCAD_EXPECTS(csr_values.size() == levels_.front().a.nonzeros());
    built_ = false;
    {
        // U (R A C) U^-1: entry (i, j) times U_i / U_j.
        Csr& fine = levels_.front().a;
        for (std::size_t i = 0; i < fine.rows; ++i) {
            for (std::size_t k = at(fine.offsets[i]); k < at(fine.offsets[i + 1]); ++k) {
                fine.values[k] = csr_values[k] * similarity_[i] / similarity_[at(fine.columns[k])];
            }
        }
    }
    levels_.resize(1);
    std::vector<double> candidate = candidate_;
    std::size_t stored = 0;
    for (std::size_t l = 0;; ++l) {
        const Csr& a = levels_[l].a;
        const bool last_allowed = l + 1 >= static_cast<std::size_t>(settings_.max_levels) ||
                                  a.rows <= settings_.coarse_size;
        if (!coarsened_) {
            if (last_allowed) break;
            const double theta_l =
                settings_.strength_threshold * std::pow(settings_.strength_decay, static_cast<double>(l));
            amg::Coarsening c = settings_.strength == StrengthMeasure::scaled && settings_.strength_decay == 1.0
                                    ? amg::coarsen(a, settings_.strength_threshold)
                                    : amg::coarsen(a, settings_.strength, theta_l, candidate);
            if (c.aggregates == 0 ||
                static_cast<double>(c.aggregates) > amg::min_coarsening * static_cast<double>(a.rows)) {
                // Stalled above coarse_size: only smoothed, unless factorize_stalled and small enough.
                coarse_smoothed_ = !(settings_.factorize_stalled &&
                                     a.rows <= amg::stalled_factor * settings_.coarse_size);
                break;
            }
            coarsening_.push_back(std::move(c));
        } else if (l >= coarsening_.size()) {
            break;
        }
        const amg::Coarsening& c = coarsening_[l];
        const std::size_t n = a.rows;
        const std::size_t nc = c.aggregates;

        // Tentative prolongator: the candidate normalized on each aggregate; the coarse candidate
        // is its norm there.
        std::vector<double> norm2(nc, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            if (c.aggregate[i] >= 0) norm2[at(c.aggregate[i])] += candidate[i] * candidate[i];
        }
        std::vector<double> coarse_candidate(nc);
        for (std::size_t j = 0; j < nc; ++j) coarse_candidate[j] = std::sqrt(norm2[j]);
        Csr p0;
        p0.rows = n;
        p0.cols = nc;
        p0.offsets.assign(n + 1, 0);
        for (std::size_t i = 0; i < n; ++i) {
            if (c.aggregate[i] >= 0) {
                const double scale = coarse_candidate[at(c.aggregate[i])];
                p0.columns.push_back(c.aggregate[i]);
                p0.values.push_back(scale > 0.0 ? candidate[i] / scale : 0.0);
            }
            p0.offsets[i + 1] = static_cast<Index>(p0.columns.size());
        }
        // Filtered matrix: strong entries and the diagonal; weak ones lumped into the diagonal
        // weighted by the candidate, so A_F b = A b for the near-null vector b.
        Csr filtered;
        filtered.rows = filtered.cols = n;
        filtered.offsets.assign(n + 1, 0);
        for (std::size_t i = 0; i < n; ++i) {
            double lumped = 0.0;
            for (std::size_t k = at(a.offsets[i]); k < at(a.offsets[i + 1]); ++k) {
                const std::size_t j = at(a.columns[k]);
                if (j == i) {
                    lumped += a.values[k];
                } else if (c.strong[k] == 0 && candidate[i] != 0.0) {
                    lumped += a.values[k] * candidate[j] / candidate[i];
                }
            }
            bool placed = false;
            for (std::size_t k = at(a.offsets[i]); k < at(a.offsets[i + 1]); ++k) {
                const std::size_t j = at(a.columns[k]);
                if (!placed && j >= i) {
                    filtered.columns.push_back(static_cast<Index>(i));
                    filtered.values.push_back(lumped);
                    placed = true;
                }
                if (j != i && c.strong[k] != 0) {
                    filtered.columns.push_back(a.columns[k]);
                    filtered.values.push_back(a.values[k]);
                }
            }
            if (!placed) {
                filtered.columns.push_back(static_cast<Index>(i));
                filtered.values.push_back(lumped);
            }
            filtered.offsets[i + 1] = static_cast<Index>(filtered.columns.size());
        }
        // P = P0 - omega D^-1 A_F P0: P0's pattern lies in A_F P0's (A_F stores every diagonal).
        const std::vector<double> inv_f = inverse_diagonal(filtered);
        const double omega = amg::prolongator_damping / spectral_radius(filtered, inv_f);
        Csr p = amg::product(filtered, p0);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = at(p.offsets[i]); k < at(p.offsets[i + 1]); ++k) {
                const double tentative =
                    (c.aggregate[i] >= 0 && p.columns[k] == c.aggregate[i]) ? p0.values[at(p0.offsets[i])]
                                                                            : 0.0;
                p.values[k] = tentative - omega * inv_f[i] * p.values[k];
            }
        }
        Csr r = amg::transpose(p);
        Csr coarse = amg::product(r, amg::product(a, p));
        stored += p.nonzeros() + r.nonzeros() + coarse.nonzeros();
        if (!all_finite(p.values) || !all_finite(coarse.values)) {
            return std::unexpected(build_error("AMG hierarchy has a non-finite entry"));
        }
        levels_[l].p = std::move(p);
        levels_[l].r = std::move(r);
        Level& next = levels_.emplace_back();
        next.a = std::move(coarse);
        candidate = std::move(coarse_candidate);
    }
    coarsened_ = true;

    // The coarsest level, by SparseLU (its pattern is fixed once the coarsening is), unless the
    // coarsening stalled above coarse_size: then that level is only smoothed.
    const Csr& last = levels_.back().a;
    if (!coarse_smoothed_ && !coarse_analyzed_) {
        std::vector<Triplet> triplets;
        triplets.reserve(last.nonzeros());
        for (std::size_t i = 0; i < last.rows; ++i) {
            for (std::size_t k = at(last.offsets[i]); k < at(last.offsets[i + 1]); ++k) {
                triplets.push_back({static_cast<Index>(i), last.columns[k], last.values[k]});
            }
        }
        auto matrix = SparseMatrix::from_triplets(static_cast<Index>(last.rows),
                                                  static_cast<Index>(last.rows), triplets);
        if (!matrix) return std::unexpected(std::move(matrix.error()));
        coarse_ = std::make_unique<SparseLuBackend<double>>();
        if (auto analyzed = coarse_->analyze(*matrix); !analyzed) {
            return std::unexpected(build_error("AMG coarsest level could not be analyzed"));
        }
        coarse_analyzed_ = true;
    }
    if (!coarse_smoothed_) {
        const auto lu = coarse_->factorize(last.values);
        if (!lu || (lu->pivots && !(lu->pivots->ratio > 0.0))) {
            return std::unexpected(build_error("AMG coarsest level could not be factorized"));
        }
        stored += lu->fill.value_or(0);
    }

    for (Level& level : levels_) {
        level.inv_d = inverse_diagonal(level.a);
        level.lambda = spectral_radius(level.a, level.inv_d);
        for (auto* v : {&level.b, &level.x, &level.residual, &level.d}) v->assign(level.a.rows, 0.0);
    }
    ++counts_.factorizations;
    built_ = true;
    return Factorization{.pivots = std::nullopt, .perturbed_pivots = std::nullopt, .fill = stored};
}

void AmgBackend::smooth(Level& level, bool from_zero) {
    const Csr& a = level.a;
    const std::size_t n = a.rows;
    const double lambda_max = level.lambda;
    const double lambda_min = lambda_max / amg::chebyshev_ratio;
    const double theta = 0.5 * (lambda_max + lambda_min);
    const double delta = 0.5 * (lambda_max - lambda_min);
    const double sigma = theta / delta;
    const auto each = [&](auto&& body) {
        base::for_chunks(pool_, n, amg::vector_chunk, [&](std::size_t begin, std::size_t end) {
            for (std::size_t i = begin; i < end; ++i) body(i);
        });
    };
    const auto residual = [&] {
        amg::multiply(pool_, a, level.x, level.residual);
        each([&](std::size_t i) { level.residual[i] = level.b[i] - level.residual[i]; });
    };
    if (from_zero) {
        each([&](std::size_t i) {
            level.d[i] = level.inv_d[i] * level.b[i] / theta;
            level.x[i] = level.d[i];
        });
    } else {
        residual();
        each([&](std::size_t i) {
            level.d[i] = level.inv_d[i] * level.residual[i] / theta;
            level.x[i] += level.d[i];
        });
    }
    double rho = 1.0 / sigma;
    for (int step = 1; step < settings_.smoother_degree; ++step) {
        residual();
        const double rho_next = 1.0 / (2.0 * sigma - rho);
        const double keep = rho_next * rho;
        const double gain = 2.0 * rho_next / delta;
        each([&](std::size_t i) {
            level.d[i] = keep * level.d[i] + gain * level.inv_d[i] * level.residual[i];
            level.x[i] += level.d[i];
        });
        rho = rho_next;
    }
}

void AmgBackend::cycle(std::size_t l) {
    Level& level = levels_[l];
    if (l + 1 == levels_.size()) {
        if (coarse_smoothed_) {
            smooth(level, true);
            smooth(level, false);
        } else {
            coarse_->solve(level.b, level.x);
        }
        return;
    }
    smooth(level, true);
    amg::multiply(pool_, level.a, level.x, level.residual);
    base::for_chunks(pool_, level.a.rows, amg::vector_chunk, [&](std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) level.residual[i] = level.b[i] - level.residual[i];
    });
    Level& next = levels_[l + 1];
    amg::multiply(pool_, level.r, level.residual, next.b);
    cycle(l + 1);
    // x += P x_c, each row in its own order.
    const Csr& p = level.p;
    base::for_chunks(pool_, p.rows, amg::row_chunk, [&](std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            double s = 0.0;
            for (std::size_t k = at(p.offsets[i]); k < at(p.offsets[i + 1]); ++k) {
                s += p.values[k] * next.x[at(p.columns[k])];
            }
            level.x[i] += s;
        }
    });
    smooth(level, false);
}

void AmgBackend::solve(std::span<const double> b, std::span<double> x) {
    NITCAD_EXPECTS(built_);
    Level& fine = levels_.front();
    NITCAD_EXPECTS(b.size() == fine.a.rows && x.size() == b.size());
    for (std::size_t i = 0; i < b.size(); ++i) fine.b[i] = similarity_[i] * b[i];
    cycle(0);
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = fine.x[i] / similarity_[i];
    ++counts_.solves;
}

SymmetrizedAmg::SymmetrizedAmg(const AmgSettings& settings, base::ThreadPool* pool)
    : pool_(pool), amg_(settings, pool) {}

std::expected<void, base::Error> SymmetrizedAmg::analyze(const SparseMatrix& a) {
    NITCAD_EXPECTS(a.rows() == a.cols() && a.rows() > 0);
    a_.rows = a_.cols = static_cast<std::size_t>(a.rows());
    a_.offsets.assign(a.row_offsets().begin(), a.row_offsets().end());
    a_.columns.assign(a.col_indices().begin(), a.col_indices().end());
    a_.values.assign(a.nonzeros(), 0.0);
    b_.assign(a_.rows, 0.0);
    x_.assign(a_.rows, 0.0);
    return amg_.analyze(a);
}

std::expected<Factorization, base::Error> SymmetrizedAmg::factorize(std::span<const double> csr_values) {
    NITCAD_EXPECTS(csr_values.size() == a_.values.size());
    std::ranges::copy(csr_values, a_.values.begin());
    symmetrization_ = amg::symmetrize(a_);
    const std::vector<double> s = amg::symmetric_part(a_, symmetrization_.d);
    amg_.set_scales(symmetrization_.d, symmetrization_.d, true);
    return amg_.factorize(s);
}

void SymmetrizedAmg::solve(std::span<const double> r, std::span<double> z) {
    NITCAD_EXPECTS(r.size() == a_.rows && z.size() == r.size());
    const std::span<const double> d = symmetrization_.d;
    base::for_chunks(pool_, r.size(), amg::vector_chunk, [&](std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) b_[i] = d[i] * r[i];
    });
    amg_.solve(b_, x_);
    base::for_chunks(pool_, z.size(), amg::vector_chunk, [&](std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) z[i] = x_[i] / d[i];
    });
}

}  // namespace NiTCAD::linalg::detail
