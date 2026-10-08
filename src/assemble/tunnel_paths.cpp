#include "NiTCAD/assemble/tunnel_paths.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::assemble {

namespace {

using Index = mesh::TensorCells::Index;

constexpr int max_steps = 100000;
constexpr int stall_steps = 64;            // steps without a new maximum of psi
constexpr double flat_field_V_cm = 1e3;    // `extend`: the end of the band bending
constexpr double extend_lengths = 10.0;    // `extend`: at most this many max_length_cm

class Tracer {
public:
    explicit Tracer(const TunnelTraceInput& in) : in_(in), c_(*in.cells), D_(c_.dimension()) {
        const std::size_t n = in.psi.size();
        gradient_.assign(n, mesh::Point{});
        for (std::size_t node = 0; node < n; ++node) {
            const Index at = c_.index(static_cast<mesh::NodeId>(node));
            for (int a = 0; a < D_; ++a) gradient_[node][a] = derivative(at, a);
        }
    }

    [[nodiscard]] std::optional<TunnelPath> trace(std::size_t start) const {
        if (!starts(start)) return std::nullopt;
        TunnelPath path;
        path.start = start;
        if (!run(start, 1.0, path.forward, path.length_cm, &path.crossing)) return std::nullopt;
        if (in_.extend) {
            std::vector<double> unused;
            (void)run(start, -1.0, path.backward, unused, nullptr);
        }
        return path;
    }

private:
    // d psi / d axis a at a node, as numpy.gradient: second-order central inside, first-order
    // one-sided at the ends.
    double derivative(const Index& at, int a) const {
        const std::span<const double> x = c_.axis(a);
        const std::size_t j = at[a], n = x.size();
        const auto f = [&](std::size_t k) {
            Index i = at;
            i[a] = k;
            return in_.psi[static_cast<std::size_t>(c_.node(i))];
        };
        if (j == 0) return (f(1) - f(0)) / (x[1] - x[0]);
        if (j == n - 1) return (f(n - 1) - f(n - 2)) / (x[n - 1] - x[n - 2]);
        const double h1 = x[j] - x[j - 1], h2 = x[j + 1] - x[j];
        return (h1 * h1 * f(j + 1) - h2 * h2 * f(j - 1) + (h2 * h2 - h1 * h1) * f(j)) /
               (h1 * h2 * (h1 + h2));
    }

    mesh::Point point(std::size_t node) const {
        const Index at = c_.index(static_cast<mesh::NodeId>(node));
        mesh::Point p{};
        for (int a = 0; a < D_; ++a) p[a] = c_.axis(a)[at[a]];
        return p;
    }

    bool traceable(const Index& cell) const {
        if (!c_.valid_cell(cell)) return false;
        for (int k = 0; k < (1 << D_); ++k) {
            Index corner = cell;
            for (int a = 0; a < D_; ++a) corner[a] += (k >> a) & 1;
            if (in_.semiconductor[static_cast<std::size_t>(c_.node(corner))] == 0) return false;
        }
        return true;
    }

    // The start criterion (header comment).
    bool starts(std::size_t i) const {
        if (in_.semiconductor[i] == 0 || in_.contact[i] != 0) return false;
        const double length = in_.max_length_cm[i];
        if (!(length > 0.0)) return false;
        double g2 = 0.0;
        for (int a = 0; a < D_; ++a) g2 += gradient_[i][a] * gradient_[i][a];
        if (g2 == 0.0) return false;
        const double R = length + c_.largest_diagonal();
        const Index at = c_.index(static_cast<mesh::NodeId>(i));
        Index lo{}, hi{};
        for (int a = 0; a < D_; ++a) {
            const std::span<const double> x = c_.axis(a);
            const double v = x[at[a]];
            lo[a] = static_cast<std::size_t>(std::lower_bound(x.begin(), x.end(), v - R) -
                                             x.begin());
            hi[a] = static_cast<std::size_t>(std::upper_bound(x.begin(), x.end(), v + R) -
                                             x.begin());
        }
        const double Ev = in_.valence[i];
        Index k = lo;
        while (true) {
            const auto j = static_cast<std::size_t>(c_.node(k));
            if (in_.semiconductor[j] != 0 && in_.end_band[j] <= Ev) {
                double d2 = 0.0;
                for (int a = 0; a < D_; ++a) {
                    const double d = c_.axis(a)[k[a]] - c_.axis(a)[at[a]];
                    d2 += d * d;
                }
                if (d2 <= R * R) return true;
            }
            int a = 0;
            for (; a < D_; ++a) {
                if (++k[a] < hi[a]) break;
                k[a] = lo[a];
            }
            if (a == D_) return false;
        }
    }

    // The multilinear interpolation of the nodal gradients at p in a cell, its components below
    // 1e-12 of the largest set to 0 (sharp): in a transverse-uniform state the transverse gradient
    // is 0 only to rounding, which would otherwise decide which neighbouring cell a path takes
    // and move it off its row (traced exactly as in 1D instead).
    mesh::Point direction(const Index& cell, const mesh::Point& p) const {
        const mesh::TensorCells::Stencil s = c_.stencil(cell, p);
        mesh::Point d{};
        for (int k = 0; k < s.count; ++k) {
            const auto node = static_cast<std::size_t>(s.nodes[static_cast<std::size_t>(k)]);
            const double w = s.weights[static_cast<std::size_t>(k)];
            for (int a = 0; a < D_; ++a) d[a] += w * gradient_[node][a];
        }
        return sharp(d);
    }

    mesh::Point sharp(mesh::Point d) const {
        double largest = 0.0;
        for (int a = 0; a < D_; ++a) largest = std::max(largest, std::abs(d[a]));
        for (int a = 0; a < D_; ++a) {
            if (std::abs(d[a]) <= 1e-12 * largest) d[a] = 0.0;
        }
        return d;
    }

    // The traceable cell at a start node along direction g: the one g points into on every axis,
    // or else the one where g, without its components out of the cell, is longest.
    std::optional<Index> first_cell(std::size_t node, const mesh::Point& g) const {
        const Index at = c_.index(static_cast<mesh::NodeId>(node));
        std::optional<Index> best;
        double best_norm = 0.0;
        for (int k = 0; k < (1 << D_); ++k) {
            Index cell = at;
            bool ok = true;
            double norm = 0.0;
            bool inside = true;
            for (int a = 0; a < D_; ++a) {
                const bool upper = ((k >> a) & 1) != 0;  // the interval above the node
                if (!upper) {
                    if (at[a] == 0) ok = false;
                    else cell[a] = at[a] - 1;
                }
                const bool out = upper ? g[a] < 0.0 : g[a] > 0.0;
                if (out) inside = false;
                else norm += g[a] * g[a];
            }
            if (!ok || !traceable(cell)) continue;
            if (inside) return cell;
            if (norm > best_norm) {
                best_norm = norm;
                best = cell;
            }
        }
        return best;
    }

    TunnelSample sample(const Index& cell, const mesh::Point& p) const {
        const mesh::TensorCells::Stencil s = c_.stencil(cell, p);
        TunnelSample t;
        t.position = p;
        t.count = s.count;
        for (int k = 0; k < s.count; ++k) {
            const auto i = static_cast<std::size_t>(k);
            t.nodes[i] = static_cast<std::size_t>(s.nodes[i]);
            t.weights[i] = s.weights[i];
        }
        return t;
    }

    double value(std::span<const double> v, const TunnelSample& s) const {
        double sum = 0.0;
        for (int k = 0; k < s.count; ++k) {
            const auto i = static_cast<std::size_t>(k);
            sum += s.weights[i] * v[s.nodes[i]];
        }
        return sum;
    }

    // Traces from `start` along sign * grad(psi). Forward (crossing set): false if the path does
    // not cross within max_length_cm.
    bool run(std::size_t start, double sign, std::vector<TunnelSample>& samples,
             std::vector<double>& lengths, std::size_t* crossing) const {
        const double l_max = in_.max_length_cm[start];
        const double Ev = in_.valence[start];
        const double margin = 0.5 * (in_.end_band[start] - Ev);
        mesh::Point p = point(start);
        mesh::Point g = sharp(gradient_[start]);
        for (int a = 0; a < D_; ++a) g[a] *= sign;
        const std::optional<Index> first = first_cell(start, g);
        if (!first) return false;
        Index cell = *first;
        samples.push_back(sample(cell, p));
        double length = 0.0, best = sign * in_.psi[start];
        int stall = 0;
        bool crossed = false;
        for (int step = 0; step < max_steps; ++step) {
            mesh::Point d = direction(cell, p);
            for (int a = 0; a < D_; ++a) d[a] *= sign;
            // Leave through the faces p sits on where d points out; a face to a cell that is not
            // traced is a Neumann boundary.
            for (int a = 0; a < D_; ++a) {
                const std::span<const double> x = c_.axis(a);
                if (d[a] > 0.0 && p[a] == x[cell[a] + 1]) {
                    Index next = cell;
                    ++next[a];
                    if (traceable(next)) cell = next;
                    else d[a] = 0.0;
                } else if (d[a] < 0.0 && p[a] == x[cell[a]]) {
                    Index next = cell;
                    if (cell[a] > 0) --next[a];
                    if (cell[a] > 0 && traceable(next)) cell = next;
                    else d[a] = 0.0;
                }
            }
            double norm = 0.0;
            for (int a = 0; a < D_; ++a) norm += d[a] * d[a];
            norm = std::sqrt(norm);
            if (norm == 0.0) break;
            if (in_.extend && (crossed || crossing == nullptr) &&
                norm * in_.V_T < flat_field_V_cm) {
                break;
            }
            // Straight to the nearer of the next face and half the cell's smallest width.
            double t_face = std::numeric_limits<double>::infinity();
            double t_axis[3] = {t_face, t_face, t_face};
            for (int a = 0; a < D_; ++a) {
                const double u = d[a] / norm;
                if (u == 0.0) continue;
                const std::span<const double> x = c_.axis(a);
                const double bound = u > 0.0 ? x[cell[a] + 1] : x[cell[a]];
                t_axis[a] = std::max((bound - p[a]) / u, 0.0);
                t_face = std::min(t_face, t_axis[a]);
            }
            const double h = 0.5 * c_.smallest_width(cell);
            mesh::Point q = p;
            // A step that would stop within 1e-6 of a step short of the face lands on it: else
            // the rounding of the potential decides between a sample just short of the face plus
            // a near-empty step onto it, and the face alone, and the geometry (its sample count
            // and crossing index) changes with the rounding.
            if (t_face <= h * (1.0 + 1e-6)) {
                for (int a = 0; a < D_; ++a) {
                    const double u = d[a] / norm;
                    if (u == 0.0) continue;
                    const std::span<const double> x = c_.axis(a);
                    if (t_axis[a] <= t_face * (1.0 + 1e-12)) {
                        q[a] = u > 0.0 ? x[cell[a] + 1] : x[cell[a]];  // land on the face
                    } else {
                        q[a] = p[a] + t_face * u;
                    }
                }
            } else {
                for (int a = 0; a < D_; ++a) {
                    const std::span<const double> x = c_.axis(a);
                    q[a] = std::clamp(p[a] + h * d[a] / norm, x[cell[a]], x[cell[a] + 1]);
                }
            }
            double L = 0.0;
            for (int a = 0; a < D_; ++a) L += (q[a] - p[a]) * (q[a] - p[a]);
            L = std::sqrt(L);
            if (L == 0.0) break;
            TunnelSample s = sample(cell, q);
            bool on_contact = false;
            for (int k = 0; k < s.count; ++k) {
                const auto i = static_cast<std::size_t>(k);
                if (s.weights[i] > 0.0 && in_.contact[s.nodes[i]] != 0) on_contact = true;
            }
            if (on_contact) break;
            samples.push_back(s);
            lengths.push_back(L);
            length += L;
            p = q;
            if (crossing != nullptr) {
                const double delta = Ev - value(in_.end_band, s);
                if (!crossed && delta >= 0.0) {
                    crossed = true;
                    *crossing = samples.size() - 2;
                }
                if (!crossed && length > l_max) return false;
                if (crossed && !in_.extend && delta >= margin) break;
            }
            if (in_.extend && length > extend_lengths * l_max) break;
            const double v = sign * value(in_.psi, s);
            if (v > best) {
                best = v;
                stall = 0;
            } else if (++stall >= stall_steps) {
                break;
            }
        }
        return crossing == nullptr || crossed;
    }

    const TunnelTraceInput& in_;
    const mesh::TensorCells& c_;
    int D_;
    std::vector<mesh::Point> gradient_;  // d psi / d axis, per node [1/cm]
};

}  // namespace

bool same_geometry(const TunnelPath& a, const TunnelPath& b, double tolerance_cm) {
    if (a.start != b.start) return false;
    const std::size_t later = std::max(a.crossing, b.crossing);
    if (later - std::min(a.crossing, b.crossing) > 1) return false;
    const std::size_t compared = later + 2;  // samples up to the later segment's end
    if (a.forward.size() < compared || b.forward.size() < compared) return false;
    // The nodes a sample reads (weight not 0), in order.
    const auto read = [](const TunnelSample& s) {
        std::vector<std::size_t> v;
        for (int k = 0; k < s.count; ++k) {
            const auto i = static_cast<std::size_t>(k);
            if (s.weights[i] != 0.0) v.push_back(s.nodes[i]);
        }
        std::ranges::sort(v);
        return v;
    };
    for (std::size_t s = 0; s < compared; ++s) {
        const TunnelSample& x = a.forward[s];
        const TunnelSample& y = b.forward[s];
        if (read(x) != read(y)) return false;
        for (int d = 0; d < 3; ++d) {
            if (std::abs(x.position[d] - y.position[d]) > tolerance_cm) return false;
        }
    }
    return true;
}

TunnelPaths trace_tunnel_paths(const TunnelTraceInput& input) {
    NITCAD_EXPECTS(input.cells != nullptr && input.V_T > 0.0);
    const std::size_t n = input.psi.size();
    NITCAD_EXPECTS(input.valence.size() == n && input.end_band.size() == n &&
                   input.semiconductor.size() == n && input.contact.size() == n &&
                   input.max_length_cm.size() == n);
    const Tracer tracer(input);
    TunnelPaths paths;
    for (std::size_t i = 0; i < n; ++i) {
        if (auto p = tracer.trace(i)) paths.paths.push_back(std::move(*p));
    }
    return paths;
}

}  // namespace NiTCAD::assemble
