// The provenance of a run (Unit 42; DECISIONS.md PV1-PV3, gate PG5): every public run records the
// build, the platform, its start and its wall time; the provenance is true, and it is not part of
// the input identity.
#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <regex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/mkl_runtime.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/results/run.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/provenance.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/solve/trace.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;

namespace {

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

// The commit the source tree is at, as git reports it now; "unknown" when the build found no git.
std::string head_revision() {
    const std::string git = NITCAD_GIT;
    if (git.empty()) return "unknown";
    // cmd /c strips the outer quotes, leaving the quoted program and directory.
    const std::string command =
        "\"\"" + git + "\" -C \"" + std::string(NITCAD_SOURCE_DIR) + "\" rev-parse HEAD\"";
    std::FILE* pipe = _popen(command.c_str(), "r");
    if (pipe == nullptr) return {};
    std::string out;
    char buffer[128];
    while (std::fgets(buffer, sizeof buffer, pipe) != nullptr) out += buffer;
    _pclose(pipe);
    while (!out.empty() && (out.back() == '\n' || out.back() == '\r')) out.pop_back();
    return out;
}

// What every run records about the build and the platform, whatever the run.
void check_provenance(const results::Provenance& p, double outer_s) {
    CHECK(p.version == "0.0.0");
    CHECK(p.revision == head_revision());
#ifdef NDEBUG
    CHECK(p.build_type == "Release");
#else
    CHECK(p.build_type == "Debug");
#endif
    CHECK(p.compiler.starts_with("MSVC 19."));
    CHECK(p.eigen_version == "5.0.1-dev");
    CHECK(p.eigen_cache_bytes == std::vector<std::int64_t>{48 * 1024, 3 * 1024 * 1024,
                                                           20 * 1024 * 1024});
    const auto mkl = linalg::loaded_mkl_runtime();
    CHECK(p.mkl_version == (mkl ? mkl->version : std::string{}));
    CHECK(p.os.starts_with("Windows "));
    CHECK(!p.cpu.empty());
    CHECK(p.logical_cpus == std::thread::hardware_concurrency());
    CHECK(std::regex_match(p.started_utc,
                           std::regex(R"(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d{3}Z)")));
    CHECK(std::isfinite(p.wall_s));
    CHECK(p.wall_s > 0.0);
    CHECK(p.wall_s <= outer_s);
}

template <class Run>
auto timed(Run run) {
    const auto start = std::chrono::steady_clock::now();
    auto r = run();
    const double s =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return std::pair{std::move(r), s};
}

}  // namespace

TEST_CASE("provenance: every public run records it") {
    const device::Device d = diode();
    const std::vector<std::vector<double>> points{{0.0, 0.0}, {0.1, 0.0}};

    SECTION("sweep_bias") {
        auto [r, s] = timed([&] { return solve::sweep_bias(d, points, {}); });
        REQUIRE(r.has_value());
        check_provenance(r->run.provenance, s);
    }
    SECTION("solve_small_signal") {
        auto [r, s] = timed([&] {
            return solve::solve_small_signal(d, points, {.frequencies_Hz = {1e6}});
        });
        REQUIRE(r.has_value());
        check_provenance(r->run.provenance, s);
    }
    SECTION("solve_transient") {
        const std::vector<solve::Waveform> w{*solve::Waveform::step(0.0, 0.1, 1e-9),
                                             *solve::Waveform::step(0.0, 0.0, 1e-9)};
        solve::TransientOptions o;
        o.t_end_s = 2e-9;
        o.adaptive = false;
        o.dt_initial_s = 1e-9;
        auto [r, s] = timed([&] { return solve::solve_transient(d, w, o); });
        REQUIRE(r.has_value());
        check_provenance(r->run.provenance, s);
    }
    SECTION("trace_bias") {
        solve::TraceOptions o;
        o.start_V = {0.0, 0.0};
        o.end_V = 0.1;
        o.step_V = 0.05;
        auto [r, s] = timed([&] { return solve::trace_bias(d, o); });
        REQUIRE(r.has_value());
        check_provenance(r->run.provenance, s);
    }
}

TEST_CASE("provenance: not part of the input identity") {
    // Two records of one input made at different times: the start times differ, the identity and
    // settings do not. A record made alone has no wall time.
    const device::Device d = diode();
    const std::vector<std::vector<double>> points{{0.0, 0.0}};
    const results::RunRecord a = solve::make_run_record(d, solve::BiasOptions{}, points);
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    const results::RunRecord b = solve::make_run_record(d, solve::BiasOptions{}, points);
    REQUIRE(a.provenance.started_utc != b.provenance.started_utc);
    CHECK(a.input_identity == b.input_identity);
    CHECK(a.settings == b.settings);
    CHECK(a.provenance.wall_s == 0.0);
    CHECK(a.provenance.revision == solve::provenance().revision);
}
