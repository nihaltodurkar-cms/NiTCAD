// Unit 57 (DECISIONS.md GE1 iii, GG2): the deck comparison. Hidden ([.deck_compare]); with
// NITCAD_COMPARE_A and NITCAD_COMPARE_B set to two directories that each hold the files of the
// same `nitcad run`s of tests/decks/ (one directory per build), every file in A must have its twin
// in B and agree with it:
// - .nitsweep, .nittran, .nitac, .nitadapt: read through the io reader, the record of the run's
//   provenance (revision, build, platform, start time, wall time) cleared, and with it the
//   adaptation report's measured times (seconds and the cost model fitted to them, which are wall
//   clock like the provenance's wall time), then written again: the two rewritten files must be
//   byte-identical;
// - .nitckpt: read against the mesh of the directory's device file; the digests, biases and every
//   field compared bit for bit, the provenance skipped;
// - .csv: line by line;
// - any other file (.nitdev, .nitmesh, ...): byte for byte.
#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/io/adapt_report.hpp"
#include "NiTCAD/io/device_file.hpp"
#include "NiTCAD/io/files.hpp"
#include "fixtures.hpp"

namespace fs = std::filesystem;
using namespace NiTCAD;

namespace {

std::string environment(const char* name) {
    char* value = nullptr;
    std::size_t length = 0;
    if (_dupenv_s(&value, &length, name) != 0 || value == nullptr) return {};
    const std::unique_ptr<char, decltype(&std::free)> owned(value, &std::free);
    return value;
}

std::string bytes_of(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

std::vector<std::string> lines_of(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    std::vector<std::string> v;
    for (std::string line; std::getline(f, line);) v.push_back(line);
    return v;
}

bool same_bits(std::span<const double> a, std::span<const double> b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::bit_cast<std::uint64_t>(a[i]) != std::bit_cast<std::uint64_t>(b[i])) return false;
    }
    return true;
}

void scrub(results::AdaptReport& r) {
    r.run.provenance = {};
    r.coarse_seconds = 0.0;
    r.verification_seconds = 0.0;
    for (results::AdaptPassRow& p : r.passes) {
        p.predicted_seconds = p.predicted_seconds_nodes = 0.0;
        p.cost_a = p.cost_b = p.cost_c = p.cost_d = 0.0;
        p.cost_mesh_per_node = p.cost_transfer_per_node = p.cost_residual = 0.0;
        p.cost_alpha = p.cost_beta = p.cost_residual_nodes = 0.0;
        p.seconds_mesh = p.seconds_transfer = p.seconds_solve = p.seconds_estimate = 0.0;
        p.seconds_estimate_factorize = p.seconds_transposed = 0.0;
        for (results::AdaptCandidateRow& c : p.candidates) c.cost = 0.0;
    }
}

// The file read, scrubbed and written again under `out`; its bytes.
std::string rewritten(const fs::path& file, const fs::path& out) {
    const std::string ext = file.extension().string();
    if (ext == ".nitsweep") {
        auto s = io::load_sweep(file);
        if (!s) FAIL(file.string() << ": " << s.error().message);
        s->run.provenance = {};
        REQUIRE(io::save_sweep(out, *s).has_value());
    } else if (ext == ".nittran") {
        auto t = io::load_transient(file);
        if (!t) FAIL(file.string() << ": " << t.error().message);
        t->run.provenance = {};
        REQUIRE(io::save_transient(out, *t).has_value());
    } else if (ext == ".nitac") {
        auto a = io::load_small_signal(file);
        if (!a) FAIL(file.string() << ": " << a.error().message);
        a->run.provenance = {};
        REQUIRE(io::save_small_signal(out, *a).has_value());
    } else {
        auto r = io::load_adapt_report(file);
        if (!r) FAIL(file.string() << ": " << r.error().message);
        scrub(*r);
        REQUIRE(io::save_adapt_report(out, *r).has_value());
    }
    return bytes_of(out);
}

void compare_checkpoints(const fs::path& a, const fs::path& b) {
    fs::path device;
    for (const auto& e : fs::directory_iterator(a.parent_path())) {
        if (e.path().extension() == ".nitdev" && e.path().stem() == a.stem()) device = e.path();
    }
    REQUIRE_FALSE(device.empty());
    const auto d = io::load_device(device);
    if (!d) FAIL(device.string() << ": " << d.error().message);
    const auto ca = io::load_checkpoint(a, d->mesh());
    const auto cb = io::load_checkpoint(b, d->mesh());
    if (!ca) FAIL(a.string() << ": " << ca.error().message);
    if (!cb) FAIL(b.string() << ": " << cb.error().message);
    CHECK(ca->mesh_digest == cb->mesh_digest);
    CHECK(ca->run_identity == cb->run_identity);
    CHECK(same_bits(ca->bias_V, cb->bias_V));
    CHECK(same_bits(ca->thermal_bias_K, cb->thermal_bias_K));
    CHECK(same_bits(ca->fields.potential_V, cb->fields.potential_V));
    CHECK(same_bits(ca->fields.n_cm3, cb->fields.n_cm3));
    CHECK(same_bits(ca->fields.p_cm3, cb->fields.p_cm3));
    CHECK(same_bits(ca->fields.temperature_K, cb->fields.temperature_K));
    CHECK(same_bits(ca->fields.electron_temperature_K, cb->fields.electron_temperature_K));
    CHECK(same_bits(ca->fields.hole_temperature_K, cb->fields.hole_temperature_K));
}

}  // namespace

TEST_CASE("deck compare: two directories of deck outputs agree (Unit 57, GG2)", "[.deck_compare]") {
    const fs::path a = environment("NITCAD_COMPARE_A"), b = environment("NITCAD_COMPARE_B");
    REQUIRE_FALSE(a.empty());
    REQUIRE_FALSE(b.empty());
    const fs::path scratch = deck_test::scratch("deck_compare");
    std::size_t files = 0;
    for (const auto& e : fs::recursive_directory_iterator(a)) {
        if (!e.is_regular_file()) continue;
        const fs::path rel = fs::relative(e.path(), a);
        const std::string ext = rel.extension().string();
        if (ext == ".deck" || ext == ".out") continue;  // the inputs and the console logs
        const fs::path twin = b / rel;
        CAPTURE(rel.string());
        REQUIRE(fs::exists(twin));
        ++files;
        if (ext == ".nitsweep" || ext == ".nittran" || ext == ".nitac" || ext == ".nitadapt") {
            CHECK(rewritten(e.path(), scratch / ("a" + ext)) == rewritten(twin, scratch / ("b" + ext)));
        } else if (ext == ".nitckpt") {
            compare_checkpoints(e.path(), twin);
        } else if (ext == ".csv") {
            CHECK(lines_of(e.path()) == lines_of(twin));
        } else {
            CHECK(bytes_of(e.path()) == bytes_of(twin));
        }
    }
    std::printf("deck compare: %zu files\n", files);
    CHECK(files > 0);
}
