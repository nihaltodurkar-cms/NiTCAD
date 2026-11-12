// The intelligent decision layer (Unit 17, slice 5; DECISIONS.md Unit 17 ID2, ID3, AR5).
#include "NiTCAD/solve/adapt_policy.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <utility>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/quantum_confinement.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::solve {
namespace {

base::Error invalid(std::string message) { return {base::ErrorCode::invalid_input, std::move(message), std::nullopt}; }

double distance(const mesh::Point& a, const mesh::Point& b) {
    return std::hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

// Nodes sorted by x, for the points within a radius.
struct XIndex {
    std::vector<std::size_t> order;
    std::vector<double> xs;
    std::span<const mesh::Point> pts;

    explicit XIndex(std::span<const mesh::Point> p) : pts(p) {
        order.resize(p.size());
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::ranges::stable_sort(order, [&](std::size_t a, std::size_t b) { return p[a][0] < p[b][0]; });
        xs.reserve(p.size());
        for (const std::size_t i : order) xs.push_back(p[i][0]);
    }
    template <class F>
    void within(const mesh::Point& c, double r, F&& f) const {
        auto k = static_cast<std::size_t>(std::ranges::lower_bound(xs, c[0] - r) - xs.begin());
        for (; k < xs.size() && xs[k] <= c[0] + r; ++k) {
            if (distance(pts[order[k]], c) <= r) f(order[k]);
        }
    }
};

struct Feature {
    mesh::Point at;
    std::array<double, 3> normal;
    double span;
    double size;
};

bool is_gate_like(device::ContactKind k) {
    return k == device::ContactKind::gate || k == device::ContactKind::electrode;
}

}  // namespace

std::expected<Prediction, base::Error> predict_sizes(const device::Device& d, std::span<const std::vector<double>> points,
                                                     const BiasOptions& options) {
    const auto contacts = d.contacts();
    for (const auto& p : points) {
        if (p.size() != contacts.size()) return std::unexpected(invalid("a bias point needs one bias per contact"));
    }
    const mesh::Mesh& m = d.mesh();
    const auto pts = m.points();
    const std::size_t nn = pts.size();
    const double T = d.temperature_K();
    const double vt = base::thermal_voltage(T);
    const double q = base::q_C;
    std::vector<std::vector<std::size_t>> adj(nn);
    for (const mesh::Edge& e : m.edges()) {
        adj[static_cast<std::size_t>(e.first)].push_back(static_cast<std::size_t>(e.second));
        adj[static_cast<std::size_t>(e.second)].push_back(static_cast<std::size_t>(e.first));
    }
    const auto semi = [&](std::size_t i) { return !d.is_insulator(static_cast<mesh::NodeId>(i)); };
    const auto net = [&](std::size_t i) { return d.net_doping(static_cast<mesh::NodeId>(i)); };
    const auto eps = [&](std::size_t i) { return base::eps0_F_per_cm * d.relative_permittivity(static_cast<mesh::NodeId>(i)); };

    // The largest difference between non-gate contacts, and per gate-like contact against the rest.
    double v_junction = 0.0;
    std::vector<double> v_gate(contacts.size(), 0.0);
    for (const auto& p : points) {
        double lo = std::numeric_limits<double>::infinity(), hi = -lo;
        for (std::size_t c = 0; c < contacts.size(); ++c) {
            if (is_gate_like(contacts[c].kind)) continue;
            lo = std::min(lo, p[c]);
            hi = std::max(hi, p[c]);
        }
        if (hi >= lo) v_junction = std::max(v_junction, hi - lo);
        for (std::size_t g = 0; g < contacts.size(); ++g) {
            if (!is_gate_like(contacts[g].kind)) continue;
            bool any = false;
            for (std::size_t c = 0; c < contacts.size(); ++c) {
                if (is_gate_like(contacts[c].kind)) continue;
                any = true;
                v_gate[g] = std::max(v_gate[g], std::abs(p[g] - p[c]));
            }
            if (!any) v_gate[g] = std::max(v_gate[g], std::abs(p[g]));
        }
    }

    const XIndex index(pts);
    std::vector<Feature> features;
    Prediction out;
    // Junctions.
    for (const mesh::Edge& e : m.edges()) {
        const auto a = static_cast<std::size_t>(e.first), b = static_cast<std::size_t>(e.second);
        if (!semi(a) || !semi(b)) continue;
        const double Na = net(a), Nb = net(b);
        if (!(Na * Nb < 0.0)) continue;
        const double t = Na / (Na - Nb);
        mesh::Point at{};
        std::array<double, 3> normal{};
        for (std::size_t k = 0; k < 3; ++k) {
            at[k] = pts[a][k] + t * (pts[b][k] - pts[a][k]);
            normal[k] = (pts[b][k] - pts[a][k]) / e.length;
        }
        const double ni = physics::intrinsic_density(d.material(e.first), T);
        const double es = eps(a);
        double N1 = std::max(std::abs(Na), ni), N2 = std::max(std::abs(Nb), ni);
        const auto width = [&](double n1, double n2) {
            const double neff = n1 * n2 / (n1 + n2);
            const double vbi = std::max(0.0, vt * std::log(n1 * n2 / (ni * ni)));
            return std::sqrt(2.0 * es * (vbi + v_junction) / (q * neff));
        };
        double W = width(N1, N2);
        for (int it = 0; it < 3; ++it) {
            index.within(at, W, [&](std::size_t i) {
                if (!semi(i)) return;
                const double N = net(i);
                if (N * Na > 0.0) N1 = std::max(N1, std::abs(N));
                if (N * Nb > 0.0) N2 = std::max(N2, std::abs(N));
            });
            W = width(N1, N2);
        }
        const double neff = N1 * N2 / (N1 + N2);
        const double LD = std::sqrt(es * vt / (q * neff));
        features.push_back({at, normal, W, 0.5 * LD});
        ++out.junction_points;
    }
    // Insulator interfaces under gates.
    const bool mlda = options.models.quantum_confinement == assemble::QuantumConfinement::mlda;
    const double lambda_q = mlda ? physics::mlda_wavelength_cm(options.models.mlda_mass_n, T) : 0.0;
    // The mean direction to the node's semiconductor neighbours off the interface itself.
    std::vector<char> on_interface(nn, 0);
    const auto inward = [&](std::size_t i) {
        std::array<double, 3> n{};
        for (const std::size_t j : adj[i]) {
            if (!semi(j) || on_interface[j] != 0) continue;
            const double L = distance(pts[i], pts[j]);
            for (std::size_t k = 0; k < 3; ++k) n[k] += (pts[j][k] - pts[i][k]) / L;
        }
        const double s = std::hypot(n[0], n[1], n[2]);
        if (s > 0.0) {
            for (double& x : n) x /= s;
        }
        return n;
    };
    const auto add_interface = [&](std::size_t i, double eps_ox, double t_ox, double v) {
        if (!(v > 0.0 && t_ox > 0.0)) return;
        const double Es = eps_ox * v / (eps(i) * t_ox);
        const double lambda = std::max(vt / Es, lambda_q);
        const std::array<double, 3> n = inward(i);
        features.push_back({pts[i], n, 5.0 * lambda, lambda});
        // The depletion region under the gate (AR5): the junction's rule, L_D / 2 across the largest
        // depletion width, its surface potential at most 2 phi_F (strong inversion) or the gate bias.
        const double N = std::abs(net(i));
        const double ni = physics::intrinsic_density(d.material(static_cast<mesh::NodeId>(i)), T);
        if (N > ni) {
            const double psi = std::min(2.0 * vt * std::log(N / ni), v);
            const double Wd = std::sqrt(2.0 * eps(i) * psi / (q * N));
            const double LD = std::sqrt(eps(i) * vt / (q * N));
            features.push_back({pts[i], n, Wd, 0.5 * LD});
        }
        ++out.interface_points;
    };
    for (std::size_t g = 0; g < contacts.size(); ++g) {
        if (contacts[g].kind != device::ContactKind::gate) continue;
        for (const mesh::NodeId v : contacts[g].nodes) on_interface[static_cast<std::size_t>(v)] = 1;
    }
    std::vector<std::size_t> electrodes;
    double v_electrode = 0.0;
    for (std::size_t c = 0; c < contacts.size(); ++c) {
        if (contacts[c].kind != device::ContactKind::electrode) continue;
        v_electrode = std::max(v_electrode, v_gate[c]);
        for (const mesh::NodeId v : contacts[c].nodes) electrodes.push_back(static_cast<std::size_t>(v));
    }
    // Semiconductor nodes next to an insulator, when electrodes exist: the interface's own nodes.
    std::vector<std::size_t> insulator_of(nn, nn);
    if (!electrodes.empty()) {
        for (std::size_t i = 0; i < nn; ++i) {
            if (!semi(i)) continue;
            for (const std::size_t j : adj[i]) {
                if (!semi(j)) {
                    insulator_of[i] = j;
                    on_interface[i] = 1;
                    break;
                }
            }
        }
    }
    for (std::size_t g = 0; g < contacts.size(); ++g) {
        if (contacts[g].kind != device::ContactKind::gate) continue;
        const device::GateStack& s = contacts[g].gate;
        for (const mesh::NodeId v : contacts[g].nodes) {
            const auto i = static_cast<std::size_t>(v);
            if (semi(i)) add_interface(i, base::eps0_F_per_cm * s.oxide_relative_permittivity, s.oxide_thickness_cm, v_gate[g]);
        }
    }
    for (std::size_t i = 0; i < nn; ++i) {
        if (insulator_of[i] == nn) continue;
        double t_ox = std::numeric_limits<double>::infinity();
        for (const std::size_t k : electrodes) t_ox = std::min(t_ox, distance(pts[i], pts[k]));
        add_interface(i, eps(insulator_of[i]), t_ox, v_electrode);
    }
    // Per node, the smallest target over the features covering it.
    out.targets.size.assign(nn, std::numeric_limits<double>::infinity());
    out.targets.normal.assign(nn, std::array<double, 3>{});
    for (const Feature& f : features) {
        index.within(f.at, f.span, [&](std::size_t i) {
            if (!semi(i) || !(f.size < out.targets.size[i])) return;
            out.targets.size[i] = f.size;
            out.targets.normal[i] = f.normal;
        });
    }
    for (const double s : out.targets.size) out.predicted_nodes += std::isfinite(s) ? 1 : 0;
    return out;
}

bool edge_exceeds(const device::Device& d, const SizeTargets& t, std::size_t e) {
    const mesh::Edge& ed = d.mesh().edges()[e];
    const auto pts = d.mesh().points();
    const auto a = static_cast<std::size_t>(ed.first), b = static_cast<std::size_t>(ed.second);
    for (const std::size_t k : {a, b}) {
        const double target = t.size[k];
        if (!std::isfinite(target)) continue;
        const auto& n = t.normal[k];
        double along = ed.length;
        if (n[0] != 0.0 || n[1] != 0.0 || n[2] != 0.0) {
            along = std::abs((pts[b][0] - pts[a][0]) * n[0] + (pts[b][1] - pts[a][1]) * n[1] + (pts[b][2] - pts[a][2]) * n[2]);
        }
        if (along > target) return true;
    }
    return false;
}

std::expected<std::vector<std::size_t>, base::Error> choose_points(std::span<const IndicatorField> fields,
                                                                   std::size_t max_points,
                                                                   std::span<const char> eligible) {
    if (fields.empty()) return std::unexpected(invalid("choosing points needs a field per point"));
    if (max_points == 0) return std::unexpected(invalid("max_points must be at least 1"));
    if (!eligible.empty() && eligible.size() != fields.size()) {
        return std::unexpected(invalid("eligibility must be one flag per point"));
    }
    const std::size_t n = fields.front().excess.size();
    for (const IndicatorField& f : fields) {
        if (f.excess.size() != n) return std::unexpected(invalid("the fields differ in size"));
    }
    const auto ok = [&](std::size_t k) { return eligible.empty() || eligible[k] != 0; };
    std::size_t last = fields.size();
    for (std::size_t k = fields.size(); k-- > 0;) {
        if (ok(k)) {
            last = k;
            break;
        }
    }
    if (last == fields.size()) return std::unexpected(invalid("no point is eligible for adaptation"));
    std::vector<std::pair<double, std::size_t>> change;
    for (std::size_t k = 1; k < last; ++k) {
        if (!ok(k)) continue;
        double s = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double d = fields[k].excess[i] - fields[k - 1].excess[i];
            s += d * d;
        }
        change.emplace_back(std::sqrt(s), k);
    }
    std::ranges::stable_sort(change, [](const auto& a, const auto& b) { return a.first > b.first || (a.first == b.first && a.second < b.second); });
    std::vector<std::size_t> at{last};
    for (const auto& [d, k] : change) {
        if (at.size() >= max_points) break;
        at.push_back(k);
    }
    std::ranges::sort(at);
    return at;
}

}  // namespace NiTCAD::solve
