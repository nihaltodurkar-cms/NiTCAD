// Files (Unit 29, io/files.hpp): meshes and sweeps round-trip bit for bit; damaged, foreign and
// mismatched files are refused; writes are atomic.
#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "NiTCAD/io/files.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "record.hpp"

using namespace NiTCAD;
using base::ErrorCode;
namespace fs = std::filesystem;

namespace {

// A fresh directory per test case.
fs::path scratch(const std::string& name) {
    const fs::path dir = fs::temp_directory_path() / "nitcad_io_test" / name;
    fs::remove_all(dir);
    fs::create_directories(dir);
    return dir;
}

bool same_bits(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::bit_cast<std::uint64_t>(a[i]) != std::bit_cast<std::uint64_t>(b[i])) return false;
    }
    return true;
}

bool same_mesh(const mesh::Mesh& a, const mesh::Mesh& b) {
    if (a.dimension() != b.dimension() || a.node_count() != b.node_count()) return false;
    for (std::size_t i = 0; i < a.node_count(); ++i) {
        if (a.points()[i] != b.points()[i] || a.volumes()[i] != b.volumes()[i]) return false;
    }
    if (a.edges().size() != b.edges().size()) return false;
    for (std::size_t e = 0; e < a.edges().size(); ++e) {
        const auto& x = a.edges()[e];
        const auto& y = b.edges()[e];
        if (x.first != y.first || x.second != y.second || x.length != y.length ||
            x.coupling_area != y.coupling_area) {
            return false;
        }
    }
    if (a.boundary().size() != b.boundary().size()) return false;
    for (std::size_t p = 0; p < a.boundary().size(); ++p) {
        const auto& x = a.boundary()[p];
        const auto& y = b.boundary()[p];
        if (x.name != y.name || x.nodes != y.nodes || x.areas != y.areas) return false;
    }
    return true;
}

std::vector<double> graded(int n) {
    std::vector<double> x{0.0};
    double h = 1e-7;
    for (int i = 0; i < n; ++i) {
        x.push_back(x.back() + h);
        h *= 1.17;
    }
    return x;
}

// A sweep with every field set, including the optional ones.
results::Sweep full_sweep() {
    results::Sweep s;
    s.run.input_identity = 0x0123456789abcdefULL;
    s.run.settings = {{"newton.tol_update", 1e-8}, {"models.srh", 1.0}};
    for (int k = 0; k < 2; ++k) {
        results::BiasPoint p;
        p.bias_V = {0.1 * k, -0.0};
        p.fields = {{1.0 / 3.0, 2.0}, {1e17, 1e3}, {1e3, 1e17}, {300.0, 300.5}};
        if (k == 1) {  // Unit 52: the carrier temperatures, absent at point 0
            p.fields.electron_temperature_K = {300.0, 1234.5678};
            p.fields.hole_temperature_K = {300.0000000000001, 301.0};
        }
        p.bands = {{0.5, 0.6}, {-0.6, -0.5}, {0.0, -0.1}, {0.0, std::nan("")}};
        p.terminal_current = {1e-3, -1e-3};
        p.terminal_current_resolution = {1e-12, 1e-12};
        p.gate_charge = {0.0, 0.0};
        p.interface_trap_charge = {-1.6e-9};
        p.edge_current_n = {1e-3};
        p.edge_current_p = {2e-4};
        p.convergence = {{{1, 0.5, 1e-2, 0.25}, {2, 1e-9, 1e-12, 1.0}}, true, 2.1e-3};
        p.tunnel_paths = {{4, 1e-6, 2e5, 1e10, true}, {5, 2e-6, 1e5, 1e8, false}};
        p.path_relocations = 2;
        p.thermal_bias_K = {300.0};
        p.thermal_contact_heat = {1e-4};
        p.bias_substeps = 3;
        p.weighted_current = {3.5e-13, -3.5e-13};
        p.weighted_current_resolution = {2.4e-19, 2.4e-19};
        s.points.push_back(p);
    }
    s.stopped = base::Error{ErrorCode::non_convergence, "no convergence at point 2",
                            base::ErrorContext{.index = 2, .value = 1e-3}};
    s.unfinished = results::ConvergenceRecord{{{1, 2.0, 1.0}}, false, std::nullopt};
    return s;
}

void check_same(const results::BiasPoint& a, const results::BiasPoint& b) {
    CHECK(same_bits(a.bias_V, b.bias_V));
    CHECK(same_bits(a.fields.potential_V, b.fields.potential_V));
    CHECK(same_bits(a.fields.n_cm3, b.fields.n_cm3));
    CHECK(same_bits(a.fields.p_cm3, b.fields.p_cm3));
    CHECK(same_bits(a.fields.temperature_K, b.fields.temperature_K));
    CHECK(same_bits(a.fields.electron_temperature_K, b.fields.electron_temperature_K));
    CHECK(same_bits(a.fields.hole_temperature_K, b.fields.hole_temperature_K));
    CHECK(same_bits(a.bands.conduction_eV, b.bands.conduction_eV));
    CHECK(same_bits(a.bands.valence_eV, b.bands.valence_eV));
    CHECK(same_bits(a.bands.electron_fermi_eV, b.bands.electron_fermi_eV));
    CHECK(same_bits(a.bands.hole_fermi_eV, b.bands.hole_fermi_eV));
    CHECK(same_bits(a.terminal_current, b.terminal_current));
    CHECK(same_bits(a.terminal_current_resolution, b.terminal_current_resolution));
    CHECK(same_bits(a.gate_charge, b.gate_charge));
    CHECK(same_bits(a.interface_trap_charge, b.interface_trap_charge));
    CHECK(same_bits(a.edge_current_n, b.edge_current_n));
    CHECK(same_bits(a.edge_current_p, b.edge_current_p));
    REQUIRE(a.convergence.iterations.size() == b.convergence.iterations.size());
    for (std::size_t i = 0; i < a.convergence.iterations.size(); ++i) {
        CHECK(a.convergence.iterations[i].iteration == b.convergence.iterations[i].iteration);
        CHECK(a.convergence.iterations[i].update == b.convergence.iterations[i].update);
        CHECK(a.convergence.iterations[i].residual == b.convergence.iterations[i].residual);
        CHECK(a.convergence.iterations[i].step == b.convergence.iterations[i].step);
    }
    CHECK(a.convergence.converged == b.convergence.converged);
    CHECK(a.convergence.smallest_pivot_ratio == b.convergence.smallest_pivot_ratio);
    REQUIRE(a.tunnel_paths.size() == b.tunnel_paths.size());
    for (std::size_t i = 0; i < a.tunnel_paths.size(); ++i) {
        CHECK(a.tunnel_paths[i].start_node == b.tunnel_paths[i].start_node);
        CHECK(a.tunnel_paths[i].length_cm == b.tunnel_paths[i].length_cm);
        CHECK(a.tunnel_paths[i].field_V_per_cm == b.tunnel_paths[i].field_V_per_cm);
        CHECK(a.tunnel_paths[i].generation_cm3_s == b.tunnel_paths[i].generation_cm3_s);
        CHECK(a.tunnel_paths[i].reached == b.tunnel_paths[i].reached);
    }
    CHECK(a.path_relocations == b.path_relocations);
    CHECK(same_bits(a.thermal_bias_K, b.thermal_bias_K));
    CHECK(same_bits(a.thermal_contact_heat, b.thermal_contact_heat));
    CHECK(a.bias_substeps == b.bias_substeps);
    CHECK(same_bits(a.weighted_current, b.weighted_current));
    CHECK(same_bits(a.weighted_current_resolution, b.weighted_current_resolution));
}

void corrupt_byte(const fs::path& path, std::size_t at) {
    std::fstream f(path, std::ios::binary | std::ios::in | std::ios::out);
    f.seekg(static_cast<std::streamoff>(at));
    char c = 0;
    f.read(&c, 1);
    c = static_cast<char>(c ^ 0x5a);
    f.seekp(static_cast<std::streamoff>(at));
    f.write(&c, 1);
}

}  // namespace

TEST_CASE("files: meshes round-trip bit for bit in 1D, 2D and 3D") {
    const fs::path dir = scratch("mesh");
    const auto x = graded(30);
    const auto y = graded(7);
    const auto z = graded(5);
    for (auto made : {mesh::make_tensor_grid(x), mesh::make_tensor_grid(x, y),
                      mesh::make_tensor_grid(x, y, z)}) {
        REQUIRE(made.has_value());
        const fs::path file = dir / ("mesh" + std::to_string(made->dimension()) + ".ntc");
        REQUIRE(io::save_mesh(file, *made).has_value());
        const auto back = io::load_mesh(file);
        REQUIRE(back.has_value());
        CHECK(same_mesh(*made, *back));
        CHECK(io::mesh_digest(*made) == io::mesh_digest(*back));
    }
    // The digest sees a one-ulp change of one coordinate.
    const auto a = *mesh::make_tensor_grid(x);
    auto x2 = x;
    x2[5] = std::nextafter(x2[5], 1.0);
    CHECK(io::mesh_digest(*mesh::make_tensor_grid(x2)) != io::mesh_digest(a));
}

TEST_CASE("files: a sweep round-trips every field bit for bit") {
    const fs::path file = scratch("sweep") / "sweep.ntc";
    const results::Sweep s = full_sweep();
    REQUIRE(io::save_sweep(file, s).has_value());
    const auto back = io::load_sweep(file);
    REQUIRE(back.has_value());
    CHECK(back->run.input_identity == s.run.input_identity);
    CHECK(back->run.settings == s.run.settings);
    REQUIRE(back->points.size() == s.points.size());
    for (std::size_t k = 0; k < s.points.size(); ++k) check_same(s.points[k], back->points[k]);
    REQUIRE(back->stopped.has_value());
    CHECK(back->stopped->code == ErrorCode::non_convergence);
    CHECK(back->stopped->message == s.stopped->message);
    CHECK(back->stopped->context->index == 2);
    CHECK(back->stopped->context->value == 1e-3);
    REQUIRE(back->unfinished.has_value());
    CHECK(back->unfinished->iterations.size() == 1);
    CHECK_FALSE(back->unfinished->smallest_pivot_ratio.has_value());

    // An empty, completed sweep: no stopping error, no unfinished history.
    const fs::path empty = file.parent_path() / "empty.ntc";
    REQUIRE(io::save_sweep(empty, results::Sweep{}).has_value());
    const auto none = io::load_sweep(empty);
    REQUIRE(none.has_value());
    CHECK(none->points.empty());
    CHECK_FALSE(none->stopped.has_value());
    CHECK_FALSE(none->unfinished.has_value());
}

TEST_CASE("files: damaged, foreign and missing files are refused") {
    const fs::path dir = scratch("damaged");
    const fs::path file = dir / "sweep.ntc";
    REQUIRE(io::save_sweep(file, full_sweep()).has_value());
    const auto size = fs::file_size(file);

    // Missing.
    CHECK(io::load_sweep(dir / "absent.ntc").error().code == ErrorCode::invalid_input);
    // Another content kind.
    const auto m = *mesh::make_tensor_grid(graded(5));
    REQUIRE(io::save_mesh(dir / "mesh.ntc", m).has_value());
    CHECK(io::load_sweep(dir / "mesh.ntc").error().code == ErrorCode::invalid_input);
    CHECK(io::load_mesh(file).error().code == ErrorCode::invalid_input);
    // A byte of the payload flipped: the checksum.
    const fs::path flipped = dir / "flipped.ntc";
    fs::copy_file(file, flipped);
    corrupt_byte(flipped, 40);
    const auto e = io::load_sweep(flipped);
    REQUIRE_FALSE(e.has_value());
    CHECK(e.error().message.find("checksum") != std::string::npos);
    // Truncated.
    const fs::path cut = dir / "cut.ntc";
    fs::copy_file(file, cut);
    fs::resize_file(cut, size - 3);
    CHECK(io::load_sweep(cut).error().code == ErrorCode::invalid_input);
    // Bytes appended after a complete file.
    const fs::path longer = dir / "longer.ntc";
    fs::copy_file(file, longer);
    {
        std::ofstream f(longer, std::ios::binary | std::ios::app);
        f << "tail";
    }
    CHECK(io::load_sweep(longer).error().code == ErrorCode::invalid_input);
    // Another format version.
    const fs::path version = dir / "version.ntc";
    fs::copy_file(file, version);
    corrupt_byte(version, 8);
    CHECK(io::load_sweep(version).error().message.find("format version") != std::string::npos);
    // Not a NiTCAD file.
    const fs::path text = dir / "text.ntc";
    {
        std::ofstream f(text);
        f << "this is not a NiTCAD file, but long enough to hold a header";
    }
    CHECK(io::load_sweep(text).error().message.find("not a NiTCAD file") != std::string::npos);
}

TEST_CASE("files: a save replaces the file atomically and leaves no temporary") {
    const fs::path dir = scratch("atomic");
    const fs::path file = dir / "sweep.ntc";
    results::Sweep first = full_sweep();
    REQUIRE(io::save_sweep(file, first).has_value());
    results::Sweep second = first;
    second.run.input_identity = 42;
    REQUIRE(io::save_sweep(file, second).has_value());
    CHECK(io::load_sweep(file)->run.input_identity == 42);
    CHECK_FALSE(fs::exists(dir / "sweep.ntc.tmp"));
    // A directory that does not exist: nothing written, an error.
    const auto e = io::save_sweep(dir / "no" / "such" / "dir" / "x.ntc", first);
    REQUIRE_FALSE(e.has_value());
    CHECK(e.error().code == ErrorCode::resource_exhausted);
}

TEST_CASE("files: a checkpoint is refused on another mesh or with the wrong node count") {
    const fs::path dir = scratch("checkpoint");
    const auto m = *mesh::make_tensor_grid(graded(4));  // 5 nodes
    results::BiasPoint p;
    p.bias_V = {0.3, 0.0};
    p.fields = {{1, 2, 3, 4, 5}, {1, 1, 1, 1, 1}, {2, 2, 2, 2, 2}, {}};
    results::RunRecord run{.input_identity = 7, .settings = {}};
    REQUIRE(io::save_checkpoint(dir / "ok.ntc", m, run, p).has_value());
    const auto c = io::load_checkpoint(dir / "ok.ntc", m);
    REQUIRE(c.has_value());
    CHECK(c->bias_V == p.bias_V);
    CHECK(c->run_identity == 7);
    CHECK(c->mesh_digest == io::mesh_digest(m));
    CHECK(same_bits(c->fields.potential_V, p.fields.potential_V));
    CHECK(c->fields.temperature_K.empty());

    const auto other = *mesh::make_tensor_grid(graded(5));
    const auto wrong_mesh = io::load_checkpoint(dir / "ok.ntc", other);
    REQUIRE_FALSE(wrong_mesh.has_value());
    CHECK(wrong_mesh.error().message.find("another mesh") != std::string::npos);

    results::BiasPoint short_fields = p;
    short_fields.fields.n_cm3.pop_back();
    REQUIRE(io::save_checkpoint(dir / "short.ntc", m, run, short_fields).has_value());
    CHECK(io::load_checkpoint(dir / "short.ntc", m).error().code == ErrorCode::invalid_input);
}

TEST_CASE("files: a sweep written before Unit 30 (no step fractions) loads with full steps") {
    // The bytes a pre-Unit-30 writer produced: a sweep whose one point's convergence record has
    // iteration, update, residual and converged, but no step field (key 6). Keys as in files.cpp.
    using io::detail::Writer;
    Writer convergence;
    const std::vector<std::int64_t> iterations{1, 2};
    convergence.i64s(1, iterations);
    convergence.f64s(2, std::vector<double>{0.5, 1e-9});
    convergence.f64s(3, std::vector<double>{1e-2, 1e-12});
    convergence.u64(4, 1);
    Writer fields;
    fields.f64s(1, std::vector<double>{0.1, 0.2});
    Writer point;
    point.f64s(1, std::vector<double>{0.3, 0.0});
    point.record(2, fields);
    point.record(10, convergence);
    Writer run;
    run.u64(1, 99);
    Writer sweep;
    sweep.record(1, run);
    sweep.record(2, point);
    // The container (files.hpp): magic, version 1, kind 2 (sweep), length, payload, FNV-1a 64.
    const auto& body = sweep.bytes();
    std::vector<std::byte> file;
    const auto put = [&](std::uint64_t v, int bytes) {
        for (int b = 0; b < bytes; ++b) {
            file.push_back(static_cast<std::byte>((v >> (8 * b)) & 0xff));
        }
    };
    for (const char c : std::string("NiTCADio")) file.push_back(static_cast<std::byte>(c));
    put(1, 4);
    put(2, 4);
    put(body.size(), 8);
    file.insert(file.end(), body.begin(), body.end());
    put(io::detail::fnv1a(body), 8);
    const fs::path path = scratch("pre_unit30") / "old.ntc";
    {
        std::ofstream f(path, std::ios::binary);
        f.write(reinterpret_cast<const char*>(file.data()),
                static_cast<std::streamsize>(file.size()));
    }
    const auto back = io::load_sweep(path);
    REQUIRE(back.has_value());
    CHECK(back->run.input_identity == 99);
    REQUIRE(back->points.size() == 1);
    const auto& it = back->points[0].convergence.iterations;
    REQUIRE(it.size() == 2);
    CHECK(it[0].step == 1.0);
    CHECK(it[1].step == 1.0);
    CHECK(it[1].update == 1e-9);
    CHECK(back->points[0].bias_substeps == 0);  // also absent: the Unit 28 default
    CHECK(back->points[0].weighted_current.empty());  // and the Unit 32 fields
}
