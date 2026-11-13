// Cross-run learning, slice 1 (Unit 17b; DECISIONS.md Unit 17 ID8, AC19 part, Unit 17b LB2, LB3, LB6,
// LB11): schema 1's features and the realized-refinement labels on the legacy diode, a history
// gathered from an adaptive run, saved and loaded exactly, the manifest and the training set built
// twice bit-identical, and the refusals.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/io/learn_files.hpp"
#include "NiTCAD/mesh/locate.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/results/learn.hpp"
#include "NiTCAD/solve/adapt.hpp"
#include "NiTCAD/solve/adapt_learn.hpp"
#include "NiTCAD/solve/bias.hpp"

using namespace NiTCAD;
using Catch::Approx;

namespace {

constexpr std::size_t F = results::learn_feature_count;

// The legacy diode (adapt_loop_test.cpp's): 2 um, abrupt at 1 um, 1e17 each side.
device::Device legacy_diode(const std::vector<double>& x) {
    auto m = mesh::make_tensor_grid(x);
    REQUIRE(m.has_value());
    auto cells = mesh::TensorCells::create(*m, x);
    REQUIRE(cells.has_value());
    const std::size_t n = m->node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m->points()[i][0] < 1e-4 ? acceptors[i] : donors[i]) = 1e17;
    auto anode = m->find_boundary("x_min")->nodes;
    auto cathode = m->find_boundary("x_max")->nodes;
    auto d = device::Device::create({.mesh = std::move(*m),
                                     .temperature_K = 300.0,
                                     .regions = {{"silicon", physics::legacy_silicon()}},
                                     .node_region = std::vector<device::RegionId>(n, 0),
                                     .donors = std::move(donors),
                                     .acceptors = std::move(acceptors),
                                     .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                  {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
                                     .cells = std::move(*cells)});
    REQUIRE(d.has_value());
    return std::move(*d);
}

std::vector<double> uniform(double L, int nodes) {
    std::vector<double> a(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) a[static_cast<std::size_t>(i)] = i == nodes - 1 ? L : L * i / (nodes - 1);
    return a;
}

solve::DeviceBuilder diode_builder(int initial_nodes) {
    return [initial_nodes](const solve::MeshRefinement& r) -> std::expected<solve::AdaptableMesh, base::Error> {
        const std::vector<double> x = r.axes[0].empty() ? uniform(2e-4, initial_nodes) : r.axes[0];
        return solve::AdaptableMesh{.device = legacy_diode(x), .kind = solve::AdaptMeshKind::tensor};
    };
}

solve::BiasOptions legacy_options() {
    solve::BiasOptions o{.models = {.hole_mobility = assemble::HoleMobility::caughey_thomas}};
    o.step_control.max_cuts = 20;
    o.weighted_current = true;
    return o;
}

const std::vector<std::vector<double>> points{{-1.0, 0.0}, {0.0, 0.0}, {0.3, 0.0}, {0.6, 0.0}};

std::filesystem::path scratch(const std::string& name) {
    const auto dir = std::filesystem::temp_directory_path() / "nitcad_adapt_learn";
    std::filesystem::create_directories(dir);
    return dir / name;
}

std::string bytes_of(const std::filesystem::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

// An adaptive run's history, gathered as the runner does (LB11).
results::AdaptHistory history_of(int initial_nodes) {
    solve::HistoryRecorder rec(points, legacy_options());
    solve::AdaptOptions o;
    o.goals = {{0}};
    o.tolerance = 1e-2;
    o.max_passes = 3;
    o.on_mesh = [&](std::size_t pass, const solve::AdaptableMesh& m) { return rec.record(pass, m); };
    auto r = solve::adapt(diode_builder(initial_nodes), points, legacy_options(), o);
    REQUIRE(r.has_value());
    REQUIRE(r->mesh.has_value());
    auto h = rec.finish(*r->mesh);
    REQUIRE(h.has_value());
    h->deck_text = "the diode, " + std::to_string(initial_nodes) + " nodes";
    h->run = r->sweep.run;
    return std::move(*h);
}

}  // namespace

TEST_CASE("learn: schema 1's features of the diode's edges (LB3)", "[adapt_learn]") {
    const device::Device d = legacy_diode(uniform(2e-4, 21));
    const auto keys = solve::material_keys(d);
    REQUIRE(keys.size() == 1);
    CHECK(keys[0].rfind("semiconductor ", 0) == 0);
    auto f = solve::edge_features(d, points, legacy_options(), keys);
    REQUIRE(f.has_value());
    REQUIRE(f->rows == 20);
    REQUIRE(f->features.size() == 20 * F);
    const auto at = [&](std::size_t e, std::size_t k) { return f->features[e * F + k]; };
    for (std::size_t e = 0; e < 20; ++e) {
        CHECK(at(e, 0) == 0.0);                                        // semiconductor
        CHECK(at(e, 1) == 0.0);                                        // the one material
        CHECK(at(e, 2) == Approx(17.0).margin(1e-12));                 // log10 1e17
        CHECK(at(e, 13) == Approx(std::log10(1e-5)).epsilon(1e-12));   // 100 nm edges
        CHECK(at(e, 14) == Approx(1.0));                               // along the junction's normal
        CHECK(at(e, 5) == 0.0);                                        // no interface
        CHECK(at(e, 16) == 6.0);
        CHECK(at(e, 10) == Approx(1.0));                               // per point: |-1 - 0|
        CHECK(at(e, 11) == 0.0);                                       // no gate
        CHECK(at(e, 12) == 0.0);                                       // a terminal current
    }
    // The junction point is edge 9's midpoint (0.95 um: the doping changes sign at node 10).
    CHECK(at(9, 4) == Approx(-12.0));
    CHECK(at(0, 4) == Approx(std::log10(0.9e-4)));
    CHECK(at(9, 3) == Approx(0.0).margin(1e-12));                     // |N| equal on both sides
    CHECK(at(0, 6) == Approx(std::log10(0.5e-5)));                    // the anode's node
    CHECK(at(0, 15) == Approx(at(0, 4) - at(0, 7)));
    // An unseen material is -1 (LB8).
    const std::vector<std::string> other{"insulator 0x1p+2"};
    auto g = solve::edge_features(d, points, legacy_options(), other);
    REQUIRE(g.has_value());
    CHECK(g->features[1] == -1.0);
    // Deterministic.
    auto again = solve::edge_features(d, points, legacy_options(), keys);
    REQUIRE(again.has_value());
    CHECK(again->features == f->features);
}

TEST_CASE("learn: the realized refinement against a final mesh (LB2)", "[adapt_learn]") {
    const device::Device start = legacy_diode(uniform(2e-4, 21));
    auto f = solve::edge_features(start, points, legacy_options(), solve::material_keys(start));
    REQUIRE(f.has_value());
    // The final mesh: [0.8, 1.2] um halved, [0.9, 1.1] quartered.
    std::vector<double> x = uniform(2e-4, 21);
    for (const double v : {0.85e-4, 0.925e-4, 0.95e-4, 0.975e-4, 1.025e-4, 1.05e-4, 1.075e-4, 1.15e-4}) x.push_back(v);
    std::ranges::sort(x);
    const device::Device fin = legacy_diode(x);
    const auto locator = mesh::PointLocator::create(*fin.cells());
    const std::vector<double> l = solve::realized_levels(*f, locator, fin.mesh());
    REQUIRE(l.size() == 20);
    for (std::size_t e = 0; e < 20; ++e) {
        const double expected = (e == 9 || e == 10) ? 2.0 : (e == 8 || e == 11) ? 1.0 : 0.0;
        CHECK(l[e] == Approx(expected).margin(1e-9));
    }
    // Against itself: nothing refined.
    const auto self = mesh::PointLocator::create(*start.cells());
    for (const double v : solve::realized_levels(*f, self, start.mesh())) CHECK(v == Approx(0.0).margin(1e-12));
}

TEST_CASE("learn: a run's history, its file, the manifest and the training set are reproducible (LB6, LB11; AC19 part)",
          "[adapt_learn]") {
    const results::AdaptHistory a = history_of(21);
    REQUIRE(a.passes.size() >= 2);
    CHECK(a.schema == results::learn_schema_version);
    for (const results::HistoryPass& p : a.passes) {
        CHECK(p.features.size() == p.rows * F);
        CHECK(p.labels.size() == p.rows);
    }
    // The final pass is labelled against itself; the first has refined edges.
    for (const double v : a.passes.back().labels) CHECK(v == Approx(0.0).margin(1e-9));
    double most = 0.0;
    for (const double v : a.passes.front().labels) most = std::max(most, v);
    CHECK(most >= 1.0);

    // Saved and loaded: equal, and saved again bit-identical.
    const auto pa = scratch("a.nithist"), pa2 = scratch("a2.nithist"), pb = scratch("b.nithist");
    REQUIRE(io::save_adapt_history(pa, a).has_value());
    auto la = io::load_adapt_history(pa);
    REQUIRE(la.has_value());
    CHECK(la->passes == a.passes);
    CHECK(la->materials == a.materials);
    CHECK(la->deck_text == a.deck_text);
    REQUIRE(io::save_adapt_history(pa2, *la).has_value());
    CHECK(bytes_of(pa) == bytes_of(pa2));
    // A second device.
    REQUIRE(io::save_adapt_history(pb, history_of(31)).has_value());

    // The manifest: sorted by digest, the same twice, saved and loaded exactly.
    const std::vector<std::filesystem::path> files{pb, pa};
    auto m1 = io::make_manifest(files);
    REQUIRE(m1.has_value());
    REQUIRE(m1->entries.size() == 2);
    CHECK(m1->entries[0].digest < m1->entries[1].digest);
    const std::vector<std::filesystem::path> reversed{pa, pb};
    auto m2 = io::make_manifest(reversed);
    REQUIRE(m2.has_value());
    CHECK(*m1 == *m2);
    const auto pm = scratch("set.manifest"), pm2 = scratch("set2.manifest");
    REQUIRE(io::save_manifest(pm, *m1).has_value());
    REQUIRE(io::save_manifest(pm2, *m2).has_value());
    CHECK(bytes_of(pm) == bytes_of(pm2));
    auto lm = io::load_manifest(pm);
    REQUIRE(lm.has_value());
    CHECK(*lm == *m1);

    // The training set, twice: identical (AC19: "built twice from the same manifest is bit-identical").
    auto t1 = io::build_training_set(*lm);
    auto t2 = io::build_training_set(*lm);
    REQUIRE(t1.has_value());
    REQUIRE(t2.has_value());
    CHECK(*t1 == *t2);
    CHECK(t1->rows == t1->labels.size());
    CHECK(t1->features.size() == t1->rows * F);
    for (const double v : t1->labels) CHECK(!std::isnan(v));
    // The stride cap (LB6): at most 5 per device, taken from the start of each.
    auto capped = io::build_training_set(*lm, {}, 5);
    REQUIRE(capped.has_value());
    CHECK(capped->rows <= 10);
    CHECK(capped->rows >= 2);

    // Refusals.
    const std::vector<std::filesystem::path> dup{pa, pa2};
    auto d = io::make_manifest(dup);
    REQUIRE(!d.has_value());
    CHECK(d.error().message.find("duplicates") != std::string::npos);
    results::Manifest tampered = *m1;
    tampered.entries[0].digest ^= 1;
    auto bad = io::build_training_set(tampered);
    REQUIRE(!bad.has_value());
    CHECK(bad.error().message.find("digest") != std::string::npos);
    {
        std::ofstream f(scratch("schema2.manifest"), std::ios::binary);
        f << "nitcad-learn-manifest 1\nschema 2\n";
    }
    auto s2 = io::load_manifest(scratch("schema2.manifest"));
    REQUIRE(!s2.has_value());
    CHECK(s2.error().message.find("schema version 2") != std::string::npos);
    results::AdaptHistory other = a;
    other.schema = 2;
    CHECK(!io::save_adapt_history(scratch("x.nithist"), other).has_value());
    {
        std::ofstream f(scratch("unsorted.manifest"), std::ios::binary);
        f << "nitcad-learn-manifest 1\nschema 1\n0000000000000002 b\n0000000000000001 a\n";
    }
    CHECK(!io::load_manifest(scratch("unsorted.manifest")).has_value());
}

// --- Slice 2: the model (LB5, LB7, LB10; ID9, ID10) ----------------------------------------------

namespace {

// A synthetic training set: n rows per device, feature 13 uniform in [0, 1) by a fixed sequence,
// feature 2 another, the label a step in feature 13 plus a small slope in feature 2.
results::TrainingSet synthetic(std::size_t devices, std::size_t n) {
    results::TrainingSet t;
    t.materials = {"m"};
    for (std::size_t d = 0; d < devices; ++d) t.devices.push_back({0x100 + d, "d" + std::to_string(d)});
    for (std::size_t d = 0; d < devices; ++d) {
        for (std::size_t i = 0; i < n; ++i) {
            const double u = std::fmod(0.6180339887 * static_cast<double>(d * n + i + 1), 1.0);
            const double v = std::fmod(0.4142135624 * static_cast<double>(d * n + i + 1), 1.0);
            std::vector<double> row(F, 0.0);
            row[13] = u;
            row[2] = v;
            t.features.insert(t.features.end(), row.begin(), row.end());
            t.labels.push_back((u > 0.5 ? 2.0 : 0.0) + 0.1 * v);
            t.device.push_back(static_cast<std::uint32_t>(d));
            t.pass.push_back(static_cast<std::uint32_t>(i % 2));
            ++t.rows;
        }
    }
    return t;
}

std::vector<std::size_t> all_rows(const results::TrainingSet& t) {
    std::vector<std::size_t> r(t.rows);
    for (std::size_t i = 0; i < t.rows; ++i) r[i] = i;
    return r;
}

}  // namespace

TEST_CASE("learn: the model's parts, exactly (LB5)", "[adapt_learn]") {
    // Cut points.
    CHECK(solve::cut_points(std::vector<double>{3, 1, 2, 2}, 64) == std::vector<double>{1, 2});
    CHECK(solve::cut_points(std::vector<double>{5, 5, 5}, 64).empty());
    {
        std::vector<double> v;
        for (int i = 0; i < 100; ++i) v.push_back(i);
        const auto c = solve::cut_points(v, 4);  // ranks floor(j 100 / 5) = 20, 40, 60, 80
        CHECK(c == std::vector<double>{20, 40, 60, 80});
    }
    // Rank correlation.
    CHECK(solve::rank_correlation(std::vector<double>{1, 2, 3}, std::vector<double>{10, 20, 30}) == Approx(1.0));
    CHECK(solve::rank_correlation(std::vector<double>{1, 2, 3}, std::vector<double>{3, 2, 1}) == Approx(-1.0));
    CHECK(solve::rank_correlation(std::vector<double>{1, 1, 1}, std::vector<double>{3, 2, 1}) == 0.0);
    CHECK(solve::rank_correlation(std::vector<double>{1, 2, 2, 3}, std::vector<double>{1, 2, 3, 4}) ==
          Approx(0.9486832981).epsilon(1e-9));
    // The default split (LB7).
    CHECK(solve::default_held_out(8) == std::vector<std::size_t>{3, 7});
    CHECK(solve::default_held_out(4) == std::vector<std::size_t>{3});
    CHECK(solve::default_held_out(2) == std::vector<std::size_t>{1});
    CHECK(solve::default_held_out(1).empty());
    // Trees learn the step; ridge cannot (ID9).
    const results::TrainingSet t = synthetic(1, 400);
    const auto rows = all_rows(t);
    const solve::BoostedTrees bt = solve::fit_boosted_trees(t.features, t.labels, rows);
    REQUIRE(bt.trees.size() == 200);
    for (const auto& tree : bt.trees) CHECK(tree.size() <= 31);  // depth 4
    CHECK(bt.base == Approx(std::accumulate(t.labels.begin(), t.labels.end(), 0.0) / 400.0));
    // Away from the step (the samples between its nearest cut and 0.5 fall on the cut's side).
    double worst = 0.0;
    for (const std::size_t i : rows) {
        if (std::abs(t.features[i * F + 13] - 0.5) < 0.02) continue;
        worst = std::max(worst, std::abs(solve::predict_trees(bt.base, bt.trees, &t.features[i * F]) - t.labels[i]));
    }
    CHECK(worst < 0.25);  // against a step of 2
    CHECK(bt.importance[13] > 10.0 * bt.importance[2]);
    for (std::size_t k = 0; k < F; ++k) {
        if (k != 13 && k != 2) CHECK(bt.importance[k] == 0.0);
    }
    // The first split of the first tree is the step's.
    CHECK(bt.trees[0][0].feature == 13);
    CHECK(bt.trees[0][0].threshold == Approx(0.5).margin(0.02));
    // A tie between two equal features goes to the lower index.
    {
        results::TrainingSet twin = t;
        for (std::size_t i = 0; i < twin.rows; ++i) twin.features[i * F + 4] = twin.features[i * F + 13];
        const auto b2 = solve::fit_boosted_trees(twin.features, twin.labels, rows);
        CHECK(b2.trees[0][0].feature == 4);
    }
    // min_leaf holds.
    {
        solve::TreeOptions o;
        o.trees = 1;
        o.min_leaf = 150;
        const auto b3 = solve::fit_boosted_trees(t.features, t.labels, rows, o);
        REQUIRE(b3.trees.size() == 1);
        std::vector<std::size_t> leaf_count(b3.trees[0].size(), 0);
        for (const std::size_t i : rows) {
            std::size_t node = 0;
            const double* row = &t.features[i * F];
            while (b3.trees[0][node].feature >= 0) {
                const auto& nd = b3.trees[0][node];
                node = static_cast<std::size_t>(row[nd.feature] <= nd.threshold ? nd.left : nd.right);
            }
            ++leaf_count[node];
        }
        for (std::size_t k = 0; k < leaf_count.size(); ++k) {
            if (b3.trees[0][k].feature < 0) CHECK(leaf_count[k] >= 150);
        }
    }
    // The ridge recovers a linear law.
    {
        results::TrainingSet lin = t;
        for (std::size_t i = 0; i < lin.rows; ++i) lin.labels[i] = 3.0 + 2.0 * lin.features[i * F + 2];
        auto r = solve::fit_ridge(lin.features, lin.labels, rows, 1e-6);
        REQUIRE(r.has_value());
        for (const std::size_t i : {std::size_t{0}, std::size_t{17}, std::size_t{399}}) {
            CHECK(solve::predict_ridge(*r, &lin.features[i * F]) == Approx(lin.labels[i]).epsilon(1e-6));
        }
    }
}

TEST_CASE("learn: training by device, its metrics, its file and the version refusals (LB7, LB10, ID10; AC19 part)",
          "[adapt_learn]") {
    // Synthetic: 4 devices, the last held out by default.
    const results::TrainingSet t = synthetic(4, 300);
    const auto held = solve::default_held_out(4);
    auto m = solve::train_predictor(t, held);
    REQUIRE(m.has_value());
    CHECK(m->training.size() == 3);
    REQUIRE(m->held_out.size() == 1);
    CHECK(m->held_out[0].path == "d3");
    CHECK(m->training_samples == 900);
    CHECK(m->metrics.samples == 300);
    CHECK(m->metrics.first_pass_samples == 150);
    CHECK(m->metrics.rank_correlation > m->metrics.ridge_rank_correlation);  // trees beat ridge on a step
    CHECK(m->metrics.rank_correlation > 0.9);
    CHECK(m->metrics.recall > 0.95);
    CHECK(m->metrics.precision > 0.95);
    CHECK(m->feature_min[13] >= 0.0);
    CHECK(m->feature_max[13] < 1.0);
    CHECK(m->categories[1] == std::vector<double>{0.0});
    CHECK(m->model_class == "gradient_boosted_trees");
    CHECK(m->format == results::learn_model_format);
    CHECK(m->schema == results::learn_schema_version);
    CHECK(!m->provenance.version.empty());
    // Trained twice: equal, and its files bit-identical (AC19).
    auto again = solve::train_predictor(t, held);
    REQUIRE(again.has_value());
    CHECK(*again == *m);
    const auto p1 = scratch("m1.nitlearn"), p2 = scratch("m2.nitlearn");
    REQUIRE(io::save_learned_model(p1, *m).has_value());
    REQUIRE(io::save_learned_model(p2, *again).has_value());
    CHECK(bytes_of(p1) == bytes_of(p2));
    auto loaded = io::load_learned_model(p1);
    REQUIRE(loaded.has_value());
    CHECK(*loaded == *m);
    for (std::size_t i = 0; i < t.rows; i += 37) {
        CHECK(solve::predict_levels(*loaded, &t.features[i * F]) == solve::predict_levels(*m, &t.features[i * F]));
    }
    // Another format or schema is refused, both versions named (ID10).
    {
        results::LearnedModel other = *m;
        other.format = 2;
        REQUIRE(io::save_learned_model(scratch("f2.nitlearn"), other).has_value());
        auto r = io::load_learned_model(scratch("f2.nitlearn"));
        REQUIRE(!r.has_value());
        CHECK(r.error().code == base::ErrorCode::invalid_input);
        CHECK(r.error().message.find("format version 2, this build's is 1") != std::string::npos);
        other = *m;
        other.schema = 7;
        REQUIRE(io::save_learned_model(scratch("s7.nitlearn"), other).has_value());
        auto s = io::load_learned_model(scratch("s7.nitlearn"));
        REQUIRE(!s.has_value());
        CHECK(s.error().message.find("feature schema version 7, this build's is 1") != std::string::npos);
    }
    // Another content kind.
    REQUIRE(io::save_adapt_history(scratch("kind.nithist"), history_of(21)).has_value());
    CHECK(!io::load_learned_model(scratch("kind.nithist")).has_value());
    // The split's refusals.
    CHECK(!solve::train_predictor(t, std::vector<std::size_t>{}).has_value());
    CHECK(!solve::train_predictor(t, std::vector<std::size_t>{0, 1, 2, 3}).has_value());
    CHECK(!solve::train_predictor(t, std::vector<std::size_t>{4}).has_value());
    CHECK(!solve::train_predictor(t, std::vector<std::size_t>{1, 1}).has_value());
    // A held-out device is never trained on: changing its labels leaves the trees unchanged.
    {
        results::TrainingSet u = t;
        for (std::size_t i = 0; i < u.rows; ++i) {
            if (u.device[i] == 3) u.labels[i] = 100.0;
        }
        auto mu = solve::train_predictor(u, held);
        REQUIRE(mu.has_value());
        CHECK(mu->trees == m->trees);
        CHECK(mu->base == m->base);
        CHECK(mu->metrics.rank_correlation != m->metrics.rank_correlation);
    }
    // The diode's histories (slice 1's) train end to end.
    {
        const auto pa = scratch("t21.nithist"), pb = scratch("t31.nithist");
        REQUIRE(io::save_adapt_history(pa, history_of(21)).has_value());
        REQUIRE(io::save_adapt_history(pb, history_of(31)).has_value());
        const std::vector<std::filesystem::path> files{pa, pb};
        auto man = io::make_manifest(files);
        REQUIRE(man.has_value());
        auto set = io::build_training_set(*man);
        REQUIRE(set.has_value());
        solve::TreeOptions small;
        small.min_leaf = 2;
        auto dm = solve::train_predictor(*set, solve::default_held_out(2), small);
        REQUIRE(dm.has_value());
        CHECK(dm->held_out.size() == 1);
        CHECK(dm->metrics.samples > 0);
    }
}

// --- Slice 3: the predictor's rules (LB4, LB8) ---------------------------------------------------

TEST_CASE("learn: the training range and the first mesh's targets, exactly (LB4, LB8)", "[adapt_learn]") {
    results::LearnedModel m;
    m.feature_min.assign(F, 0.0);
    m.feature_max.assign(F, 1.0);
    m.categories.assign(F, {});
    m.categories[0] = {0.0};
    m.categories[1] = {0.0};
    std::vector<double> row(F, 0.5);
    row[0] = row[1] = 0.0;
    CHECK(!solve::out_of_range(m, row.data()));
    row[13] = 1.04;  // inside the 5% margin
    CHECK(!solve::out_of_range(m, row.data()));
    row[13] = 1.06;
    CHECK(solve::out_of_range(m, row.data()));
    row[13] = -0.04;
    CHECK(!solve::out_of_range(m, row.data()));
    row[13] = -0.06;
    CHECK(solve::out_of_range(m, row.data()));
    row[13] = 0.5;
    row[1] = 1.0;  // an unseen material
    CHECK(solve::out_of_range(m, row.data()));
    row[1] = -1.0;  // a material not in the model's list
    CHECK(solve::out_of_range(m, row.data()));
    row[1] = 0.0;
    row[4] = std::numeric_limits<double>::quiet_NaN();
    CHECK(solve::out_of_range(m, row.data()));
    // The fraction over a device's edges, against LB8's limit.
    solve::EdgeFeatures e;
    e.rows = 20;
    e.features.assign(20 * F, 0.5);
    for (std::size_t r = 0; r < 20; ++r) e.features[r * F] = e.features[r * F + 1] = 0.0;
    CHECK(solve::out_of_range_fraction(m, e) == 0.0);
    e.features[13] = 2.0;  // one edge of 20: 5%, not above the limit
    CHECK(solve::out_of_range_fraction(m, e) == Approx(0.05));
    CHECK(!(solve::out_of_range_fraction(m, e) > solve::out_of_range_limit));
    e.features[F + 13] = 2.0;
    CHECK(solve::out_of_range_fraction(m, e) > solve::out_of_range_limit);

    // Targets: a model predicting l^ = 1 everywhere halves every edge, normal along it; 0.4 asks nothing.
    const device::Device d = legacy_diode(uniform(2e-4, 21));
    results::LearnedModel one;
    one.materials = solve::material_keys(d);
    one.base = 1.0;
    auto t = solve::learned_targets(one, d, points, legacy_options());
    REQUIRE(t.has_value());
    REQUIRE(t->size.size() == 21);
    for (std::size_t i = 0; i < 21; ++i) {
        CHECK(t->size[i] == Approx(0.5e-5));
        CHECK(t->normal[i][0] == Approx(1.0));
    }
    one.base = 0.4;
    auto none = solve::learned_targets(one, d, points, legacy_options());
    REQUIRE(none.has_value());
    for (const double s : none->size) CHECK(std::isinf(s));
    one.base = 0.5;  // the threshold itself refines
    auto half = solve::learned_targets(one, d, points, legacy_options());
    REQUIRE(half.has_value());
    CHECK(half->size[0] == Approx(1e-5 * std::exp2(-0.5)));
}
