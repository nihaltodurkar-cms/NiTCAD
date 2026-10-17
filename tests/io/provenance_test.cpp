// Provenance in files (Unit 42; DECISIONS.md PV4, gates PG2 and PG3): it round-trips in sweeps and
// checkpoints, and files written before Unit 42 keep loading, with an empty provenance.
#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <string>
#include <vector>

#include "NiTCAD/io/device_file.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;

namespace {

fs::path scratch(const std::string& name) {
    const fs::path dir = fs::temp_directory_path() / "nitcad_io_test" / name;
    fs::remove_all(dir);
    fs::create_directories(dir);
    return dir;
}

// Every field set, none to its default.
results::Provenance full() {
    return {.version = "1.2.3",
            .revision = "0123456789abcdef0123456789abcdef01234567",
            .dirty = true,
            .build_type = "Release",
            .compiler = "MSVC 19.51.36260.0",
            .eigen_version = "5.0.1-dev",
            .eigen_cache_bytes = {49152, 3145728, 20971520},
            .mkl_version = "Intel(R) oneAPI Math Kernel Library Version 2025.2",
            .os = "Windows 10.0.26300",
            .cpu = "Intel(R) Core(TM) Ultra 5 225F",
            .logical_cpus = 10,
            .started_utc = "2026-10-17T09:12:03.123Z",
            .wall_s = 12.5};
}

const fs::path pre_unit42 = fs::path(NITCAD_IO_FIXTURES) / "pre_unit42";

}  // namespace

TEST_CASE("provenance: round-trips in sweeps and checkpoints") {
    const fs::path dir = scratch("provenance_round_trip");
    results::Sweep sweep;
    sweep.run = {.input_identity = 7, .settings = {{"newton.tol_update", 1e-8}},
                 .provenance = full()};
    REQUIRE(io::save_sweep(dir / "s.nitsweep", sweep).has_value());
    const auto s = io::load_sweep(dir / "s.nitsweep");
    REQUIRE(s.has_value());
    CHECK(s->run.provenance == full());
    CHECK(s->run.input_identity == 7);
    CHECK(s->run.settings == sweep.run.settings);

    // The empty provenance (a record made by hand) round-trips as empty.
    sweep.run.provenance = {};
    REQUIRE(io::save_sweep(dir / "e.nitsweep", sweep).has_value());
    const auto e = io::load_sweep(dir / "e.nitsweep");
    REQUIRE(e.has_value());
    CHECK(e->run.provenance == results::Provenance{});

    const std::vector<double> x{0.0, 1e-4, 2e-4};
    const mesh::Mesh m = *mesh::make_tensor_grid(x);
    results::BiasPoint point;
    point.bias_V = {0.1, 0.0};
    point.fields = {.potential_V = {0.0, 0.1, 0.2}, .n_cm3 = {1, 2, 3}, .p_cm3 = {4, 5, 6}};
    const results::RunRecord run{.input_identity = 9, .settings = {}, .provenance = full()};
    REQUIRE(io::save_checkpoint(dir / "c.nitckpt", m, run, point).has_value());
    const auto c = io::load_checkpoint(dir / "c.nitckpt", m);
    REQUIRE(c.has_value());
    CHECK(c->provenance == full());
    CHECK(c->run_identity == 9);
}

TEST_CASE("provenance: files written before Unit 42 load with none") {
    // Written by nitcad built from main at 67e8fe5 (bar.deck beside them).
    const auto sweep = io::load_sweep(pre_unit42 / "bar.nitsweep");
    REQUIRE(sweep.has_value());
    CHECK(sweep->run.provenance == results::Provenance{});
    CHECK(sweep->run.input_identity != 0);
    CHECK(!sweep->run.settings.empty());
    REQUIRE(sweep->points.size() == 2);
    CHECK(sweep->points[1].bias_V == std::vector<double>{0.01, 0.0});
    CHECK(sweep->points[1].fields.potential_V.size() == 11);

    const auto mesh = io::load_mesh(pre_unit42 / "bar.nitmesh");
    REQUIRE(mesh.has_value());
    const auto c = io::load_checkpoint(pre_unit42 / "bar.nitckpt", *mesh);
    REQUIRE(c.has_value());
    CHECK(c->provenance == results::Provenance{});
    CHECK(c->run_identity == sweep->run.input_identity);
    CHECK(c->bias_V == std::vector<double>{0.01, 0.0});
    CHECK(c->fields.potential_V == sweep->points[1].fields.potential_V);

    const auto device = io::load_device(pre_unit42 / "bar.nitdev");
    REQUIRE(device.has_value());
    CHECK(device->mesh().node_count() == 11);

    // Loaded and saved again, a pre-Unit-42 sweep keeps every value (and gains an empty
    // provenance record).
    const fs::path dir = scratch("provenance_resave");
    REQUIRE(io::save_sweep(dir / "again.nitsweep", *sweep).has_value());
    const auto again = io::load_sweep(dir / "again.nitsweep");
    REQUIRE(again.has_value());
    CHECK(again->run.settings == sweep->run.settings);
    CHECK(again->points[1].fields.potential_V == sweep->points[1].fields.potential_V);
    CHECK(again->points[1].terminal_current == sweep->points[1].terminal_current);
}
