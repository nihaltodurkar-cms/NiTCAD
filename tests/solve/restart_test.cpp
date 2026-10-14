// Checkpoint and restart (Unit 29): a sweep stopped part-way and restarted from a checkpoint file
// reaches the same points as the uninterrupted sweep; a saved sweep reads back unchanged.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;

namespace {

// The Unit 10 sweep fixture: the legacy 1e17/1e17 silicon diode.
device::Device diode() {
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        (m.points()[i][0] < 1e-4 ? acceptors : donors)[i] = 1e17;
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
}

double relative(const std::vector<double>& u, const std::vector<double>& v) {
    REQUIRE(u.size() == v.size());
    double worst = 0.0;
    for (std::size_t i = 0; i < u.size(); ++i) {
        worst = std::max(worst, std::abs(u[i] - v[i]) /
                                    std::max({std::abs(u[i]), std::abs(v[i]), 1e-300}));
    }
    return worst;
}

fs::path scratch() {
    const fs::path dir = fs::temp_directory_path() / "nitcad_restart_test";
    fs::remove_all(dir);
    fs::create_directories(dir);
    return dir;
}

}  // namespace

TEST_CASE("restart: a sweep resumed from a checkpoint file reaches the uninterrupted points") {
    const auto d = diode();
    std::vector<std::vector<double>> all;
    for (int k = 0; k <= 6; ++k) all.push_back({0.1 * k, 0.0});
    const auto full = solve::sweep_bias(d, all);
    REQUIRE(full.has_value());
    REQUIRE_FALSE(full->stopped.has_value());

    // The first four points, then a checkpoint of the last.
    const std::vector<std::vector<double>> head(all.begin(), all.begin() + 4);
    const auto first = solve::sweep_bias(d, head);
    REQUIRE(first.has_value());
    const fs::path file = scratch() / "at_0.3V.ntc";
    REQUIRE(io::save_checkpoint(file, d.mesh(), first->run, first->points.back()).has_value());

    // A fresh start: load it and continue.
    const auto c = io::load_checkpoint(file, d.mesh());
    REQUIRE(c.has_value());
    CHECK(c->bias_V == head.back());
    CHECK(c->run_identity == first->run.input_identity);
    const std::vector<std::vector<double>> tail(all.begin() + 4, all.end());
    solve::BiasOptions options;
    options.step_control.initial_bias_V = c->bias_V;  // read only with the step control on
    const auto rest = solve::sweep_bias(d, tail, options, &c->fields);
    REQUIRE(rest.has_value());
    REQUIRE_FALSE(rest->stopped.has_value());
    REQUIRE(rest->points.size() == tail.size());
    for (std::size_t k = 0; k < tail.size(); ++k) {
        const auto& a = rest->points[k];
        const auto& b = full->points[4 + k];
        REQUIRE(a.bias_V == b.bias_V);
        // The uninterrupted sweep carries its state in scaled units, the restart through volts and
        // cm^-3 (as the Unit 10 chain of one-point solves), so they agree to rounding of that.
        CHECK(relative(a.fields.potential_V, b.fields.potential_V) < 1e-9);
        CHECK(relative(a.fields.n_cm3, b.fields.n_cm3) < 1e-9);
        CHECK(relative(a.fields.p_cm3, b.fields.p_cm3) < 1e-9);
        CHECK(relative(a.terminal_current, b.terminal_current) < 1e-9);
    }
}

TEST_CASE("restart: a solved sweep saved and loaded is unchanged, and a loaded checkpoint state "
          "restarts exactly like the in-memory one") {
    const auto d = diode();
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 3; ++k) points.push_back({0.1 * k, 0.0});
    const auto s = solve::sweep_bias(d, points);
    REQUIRE(s.has_value());
    const fs::path dir = scratch();
    REQUIRE(io::save_sweep(dir / "sweep.ntc", *s).has_value());
    const auto back = io::load_sweep(dir / "sweep.ntc");
    REQUIRE(back.has_value());
    REQUIRE(back->points.size() == s->points.size());
    for (std::size_t k = 0; k < points.size(); ++k) {
        CHECK(back->points[k].fields.potential_V == s->points[k].fields.potential_V);
        CHECK(back->points[k].terminal_current == s->points[k].terminal_current);
        CHECK(back->points[k].edge_current_n == s->points[k].edge_current_n);
    }
    CHECK(back->run.input_identity == s->run.input_identity);

    // Restarting from the file's state and from the in-memory state is the same computation.
    REQUIRE(io::save_checkpoint(dir / "c.ntc", d.mesh(), s->run, s->points.back()).has_value());
    const auto c = io::load_checkpoint(dir / "c.ntc", d.mesh());
    REQUIRE(c.has_value());
    const std::vector<std::vector<double>> next{{0.4, 0.0}};
    const auto from_file = solve::sweep_bias(d, next, {}, &c->fields);
    const auto from_memory = solve::sweep_bias(d, next, {}, &s->points.back().fields);
    REQUIRE(from_file.has_value());
    REQUIRE(from_memory.has_value());
    CHECK(from_file->points.front().fields.potential_V ==
          from_memory->points.front().fields.potential_V);
    CHECK(from_file->points.front().terminal_current ==
          from_memory->points.front().terminal_current);
    CHECK(from_file->run.input_identity == from_memory->run.input_identity);
}
