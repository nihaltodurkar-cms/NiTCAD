// Cross-run learning: features, labels, history (Unit 17b, slice 1; DECISIONS.md LB2, LB3, LB11).
#include "NiTCAD/solve/adapt_learn.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <format>
#include <limits>
#include <numeric>
#include <optional>
#include <utility>
#include <variant>

#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/adapt_policy.hpp"

namespace NiTCAD::solve {
namespace {

base::Error invalid(std::string message) { return {base::ErrorCode::invalid_input, std::move(message), std::nullopt}; }

constexpr std::size_t F = results::learn_feature_count;

struct Fnv {
    std::uint64_t h = 0xcbf29ce484222325ull;
    void bytes(const void* p, std::size_t n) {
        const auto* b = static_cast<const unsigned char*>(p);
        for (std::size_t i = 0; i < n; ++i) {
            h ^= b[i];
            h *= 0x100000001b3ull;
        }
    }
    void u64(std::uint64_t v) {
        unsigned char b[8];
        for (int i = 0; i < 8; ++i) b[i] = static_cast<unsigned char>(v >> (8 * i));
        bytes(b, 8);
    }
    void f64(double v) { u64(std::bit_cast<std::uint64_t>(v)); }
};

double distance(const mesh::Point& a, const mesh::Point& b) { return std::hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]); }

// The nearest of a point set (lowest index on ties): points sorted by x, scanned outward from the
// query's x until the x gap alone exceeds the best distance.
class Nearest {
public:
    explicit Nearest(std::vector<mesh::Point> pts) : pts_(std::move(pts)) {
        order_.resize(pts_.size());
        std::iota(order_.begin(), order_.end(), std::size_t{0});
        std::ranges::stable_sort(order_, [&](std::size_t a, std::size_t b) { return pts_[a][0] < pts_[b][0]; });
        xs_.reserve(order_.size());
        for (const std::size_t i : order_) xs_.push_back(pts_[i][0]);
    }
    [[nodiscard]] bool empty() const noexcept { return pts_.empty(); }
    // The index of the nearest point and its distance; empty: (npos, infinity).
    [[nodiscard]] std::pair<std::size_t, double> query(const mesh::Point& p) const {
        std::size_t best = static_cast<std::size_t>(-1);
        double bd = std::numeric_limits<double>::infinity();
        const auto consider = [&](std::size_t k) {
            const std::size_t i = order_[k];
            const double d = distance(pts_[i], p);
            if (d < bd || (d == bd && i < best)) {
                bd = d;
                best = i;
            }
        };
        const auto start = static_cast<std::size_t>(std::ranges::lower_bound(xs_, p[0]) - xs_.begin());
        for (std::size_t k = start; k < xs_.size() && xs_[k] - p[0] <= bd; ++k) consider(k);
        for (std::size_t k = start; k-- > 0 && p[0] - xs_[k] <= bd;) consider(k);
        return {best, bd};
    }

private:
    std::vector<mesh::Point> pts_;
    std::vector<std::size_t> order_;
    std::vector<double> xs_;
};

}  // namespace

std::string material_key(const device::Region& region) {
    if (const auto* ins = std::get_if<physics::Insulator>(&region.material)) {
        return std::format("insulator {:a}", ins->parameters().eps_r);
    }
    const auto& p = std::get<physics::Semiconductor>(region.material).parameters();
    return std::format("semiconductor {:a} {:a} {:a} {:a}", p.eps_r, p.Eg0_eV, p.Nc300, p.Nv300);
}

std::vector<std::string> material_keys(const device::Device& d) {
    std::vector<std::string> keys;
    for (const device::Region& r : d.regions()) keys.push_back(material_key(r));
    std::ranges::sort(keys);
    keys.erase(std::ranges::unique(keys).begin(), keys.end());
    return keys;
}

std::uint64_t device_digest(const device::Device& d) {
    Fnv h;
    h.u64(mesh::graph_digest(d.mesh()));
    for (const device::RegionId r : d.node_region()) h.u64(static_cast<std::uint64_t>(static_cast<std::int64_t>(r)));
    for (const double v : d.donors()) h.f64(v);
    for (const double v : d.acceptors()) h.f64(v);
    for (const device::Contact& c : d.contacts()) {
        h.bytes(c.name.data(), c.name.size());
        h.u64(static_cast<std::uint64_t>(c.kind));
        for (const mesh::NodeId v : c.nodes) h.u64(static_cast<std::uint64_t>(static_cast<std::int64_t>(v)));
    }
    return h.h;
}

std::expected<EdgeFeatures, base::Error> edge_features(const device::Device& d, std::span<const std::vector<double>> points,
                                                       const BiasOptions& options, std::span<const std::string> materials) {
    auto fp = feature_points(d, points, options);
    if (!fp) return std::unexpected(fp.error());
    const mesh::Mesh& m = d.mesh();
    const auto pts = m.points();
    const auto& edges = m.edges();
    const double T = d.temperature_K();
    const auto regions = d.regions();
    const auto node_region = d.node_region();

    // Per region: its material's index in `materials` (-1 absent).
    std::vector<double> region_material(regions.size(), -1.0);
    for (std::size_t r = 0; r < regions.size(); ++r) {
        const std::string key = material_key(regions[r]);
        const auto it = std::ranges::lower_bound(materials, key);
        if (it != materials.end() && *it == key) region_material[r] = static_cast<double>(it - materials.begin());
    }
    const auto semi = [&](std::size_t i) { return !d.is_insulator(static_cast<mesh::NodeId>(i)); };
    // log10 max(|N|, n_i) at a semiconductor node.
    const auto log_doping = [&](std::size_t i) {
        const auto v = static_cast<mesh::NodeId>(i);
        const double ni = physics::intrinsic_density(d.material(v), T);
        return std::log10(std::max(std::abs(d.net_doping(v)), ni));
    };

    std::vector<mesh::Point> junction_at, interface_at;
    std::vector<const FeaturePoint*> junctions, interfaces;
    for (const FeaturePoint& f : fp->points) {
        if (f.kind == FeaturePoint::Kind::junction) {
            junctions.push_back(&f);
            junction_at.push_back(f.at);
        } else if (f.kind == FeaturePoint::Kind::interface) {
            interfaces.push_back(&f);
            interface_at.push_back(f.at);
        }
    }
    std::vector<mesh::Point> contact_at;
    {
        std::vector<mesh::NodeId> nodes;
        for (const device::Contact& c : d.contacts()) nodes.insert(nodes.end(), c.nodes.begin(), c.nodes.end());
        std::ranges::sort(nodes);
        nodes.erase(std::ranges::unique(nodes).begin(), nodes.end());
        for (const mesh::NodeId v : nodes) contact_at.push_back(pts[static_cast<std::size_t>(v)]);
    }
    const Nearest near_junction(std::move(junction_at)), near_interface(std::move(interface_at)),
        near_contact(std::move(contact_at));

    EdgeFeatures out;
    out.rows = edges.size();
    out.features.assign(out.rows * F, 0.0);
    out.midpoint.resize(out.rows);
    out.direction.resize(out.rows);
    out.length.resize(out.rows);
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto a = static_cast<std::size_t>(edges[e].first), b = static_cast<std::size_t>(edges[e].second);
        const double L = edges[e].length;
        mesh::Point mid{};
        std::array<double, 3> dir{};
        for (std::size_t k = 0; k < 3; ++k) {
            mid[k] = 0.5 * (pts[a][k] + pts[b][k]);
            dir[k] = (pts[b][k] - pts[a][k]) / L;
        }
        out.midpoint[e] = mid;
        out.direction[e] = dir;
        out.length[e] = L;
        double* f = &out.features[e * F];
        const bool sa = semi(a), sb = semi(b);
        f[0] = (!sa && !sb) ? 1.0 : 0.0;
        const std::size_t owner = (sa || !sb) ? a : b;
        f[1] = region_material[static_cast<std::size_t>(node_region[owner])];
        if (sa && sb) {
            const double la = log_doping(a), lb = log_doping(b);
            f[2] = 0.5 * (la + lb);
            f[3] = std::abs(lb - la) / (L * 1e4);
        } else if (sa || sb) {
            f[2] = log_doping(sa ? a : b);
        }
        const auto [jk, jd] = near_junction.query(mid);
        const auto [ik, id] = near_interface.query(mid);
        const auto [ck, cd] = near_contact.query(mid);
        static_cast<void>(ck);
        const FeaturePoint* jf = near_junction.empty() ? nullptr : junctions[jk];
        const FeaturePoint* inf = near_interface.empty() ? nullptr : interfaces[ik];
        // A distance of 0 (a midpoint on a feature point) is held at 1e-12 cm for its logarithm.
        const auto lg = [](double x) { return std::log10(std::max(x, 1e-12)); };
        f[4] = jf ? lg(jd) : 0.0;
        f[5] = inf ? lg(id) : 0.0;
        f[6] = near_contact.empty() ? 0.0 : lg(cd);
        f[7] = jf ? lg(jf->span) : 0.0;
        f[8] = jf ? lg(2.0 * jf->size) : 0.0;
        f[9] = inf ? lg(inf->size) : 0.0;
        f[10] = fp->v_junction;
        f[11] = fp->v_gate;
        f[12] = 0.0;
        f[13] = lg(L);
        const FeaturePoint* closest = jf;
        if (inf && (!jf || id < jd)) closest = inf;
        if (closest) {
            const auto& n = closest->normal;
            f[14] = std::abs(dir[0] * n[0] + dir[1] * n[1] + dir[2] * n[2]);
        }
        f[15] = jf ? lg(jd / jf->span) : 6.0;
        f[16] = inf ? lg(id / inf->size) : 6.0;
    }
    return out;
}

std::vector<double> realized_levels(const EdgeFeatures& edges, const mesh::PointLocator& locator,
                                    const mesh::Mesh& final_mesh) {
    const auto pts = final_mesh.points();
    std::vector<double> out(edges.rows, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t e = 0; e < edges.rows; ++e) {
        const mesh::Interpolant in = locator.locate(edges.midpoint[e]);
        if (!in.inside || in.count < 2) continue;
        const auto& u = edges.direction[e];
        double lo = std::numeric_limits<double>::infinity(), hi = -lo;
        for (int k = 0; k < in.count; ++k) {
            const auto& p = pts[static_cast<std::size_t>(in.nodes[static_cast<std::size_t>(k)])];
            const double s = p[0] * u[0] + p[1] * u[1] + p[2] * u[2];
            lo = std::min(lo, s);
            hi = std::max(hi, s);
        }
        const double h = hi - lo;
        if (!(h > 0.0)) continue;
        out[e] = std::max(0.0, std::log2(edges.length[e] / h));
    }
    return out;
}

std::expected<mesh::PointLocator, base::Error> adaptable_locator(const AdaptableMesh& m) {
    if (m.kind == AdaptMeshKind::tensor) {
        if (m.device.cells() == nullptr) return std::unexpected(invalid("a tensor mesh needs its cells"));
        return mesh::PointLocator::create(*m.device.cells());
    }
    if (!m.generated) return std::unexpected(invalid("a generated mesh needs its simplices"));
    return mesh::PointLocator::create(m.generated->simplices);
}

void stamp_predictor(results::RunRecord& run, std::uint64_t file_digest, const results::LearnedModel& model, bool used) {
    Fnv h;
    h.h = run.input_identity;
    const std::string_view tag = "adapt.predictor";
    h.bytes(tag.data(), tag.size());
    h.u64(file_digest);
    h.u64(model.format);
    h.u64(model.schema);
    h.u64(used ? 1 : 0);
    run.input_identity = h.h;
    run.settings.emplace_back("adapt.predictor.digest_high", static_cast<double>(file_digest >> 32));
    run.settings.emplace_back("adapt.predictor.digest_low", static_cast<double>(file_digest & 0xffffffffull));
    run.settings.emplace_back("adapt.predictor.format", static_cast<double>(model.format));
    run.settings.emplace_back("adapt.predictor.schema", static_cast<double>(model.schema));
    run.settings.emplace_back("adapt.predictor.used", used ? 1.0 : 0.0);
}

HistoryRecorder::HistoryRecorder(std::vector<std::vector<double>> points, BiasOptions options)
    : points_(std::move(points)), options_(std::move(options)) {}

std::expected<void, base::Error> HistoryRecorder::record(std::size_t pass, const AdaptableMesh& m) {
    if (pass != passes_.size()) return std::unexpected(invalid("history passes must be recorded in order"));
    Pass p;
    p.device_digest = device_digest(m.device);
    p.mesh_digest = mesh::graph_digest(m.device.mesh());
    p.materials = material_keys(m.device);
    auto f = edge_features(m.device, points_, options_, p.materials);
    if (!f) return std::unexpected(f.error());
    p.edges = std::move(*f);
    passes_.push_back(std::move(p));
    return {};
}

std::expected<results::AdaptHistory, base::Error> HistoryRecorder::finish(const AdaptableMesh& final_mesh) const {
    if (passes_.empty()) return std::unexpected(invalid("a history needs at least one solved pass"));
    auto locator = adaptable_locator(final_mesh);
    if (!locator) return std::unexpected(locator.error());
    results::AdaptHistory h;
    for (const Pass& p : passes_) h.materials.insert(h.materials.end(), p.materials.begin(), p.materials.end());
    std::ranges::sort(h.materials);
    h.materials.erase(std::ranges::unique(h.materials).begin(), h.materials.end());
    for (const Pass& p : passes_) {
        results::HistoryPass hp;
        hp.device_digest = p.device_digest;
        hp.mesh_digest = p.mesh_digest;
        hp.rows = p.edges.rows;
        hp.features = p.edges.features;
        for (std::size_t e = 0; e < p.edges.rows; ++e) {
            double& c = hp.features[e * F + 1];
            if (c < 0.0) continue;
            const std::string& key = p.materials[static_cast<std::size_t>(c)];
            c = static_cast<double>(std::ranges::lower_bound(h.materials, key) - h.materials.begin());
        }
        hp.labels = realized_levels(p.edges, *locator, final_mesh.device.mesh());
        h.passes.push_back(std::move(hp));
    }
    return h;
}

}  // namespace NiTCAD::solve
