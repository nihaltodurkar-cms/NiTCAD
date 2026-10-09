#include "NiTCAD/io/files.hpp"

#include <array>
#include <bit>
#include <cstddef>
#include <fstream>
#include <iterator>
#include <limits>
#include <optional>
#include <string>
#include <system_error>
#include <utility>

#include "record.hpp"

// Unwraps an expected into `var`, or returns its error from the enclosing function (local to this
// file: every read is a chain of fallible field reads).
#define NITCAD_IO_TRY(var, expr)                                                  \
    auto var##_expected = (expr);                                                 \
    if (!var##_expected) return std::unexpected(std::move(var##_expected.error())); \
    auto var = std::move(*var##_expected)

namespace NiTCAD::io {

namespace {

using detail::corrupt;
using detail::Reader;
using detail::Writer;

constexpr std::array<char, 8> magic{'N', 'i', 'T', 'C', 'A', 'D', 'i', 'o'};
constexpr std::size_t header_size = 8 + 4 + 4 + 8;

base::Error invalid(std::string message) {
    return {base::ErrorCode::invalid_input, std::move(message), std::nullopt};
}

base::Error unwritable(const std::filesystem::path& path, const std::string& what) {
    return {base::ErrorCode::resource_exhausted,
            "cannot write '" + path.string() + "': " + what, std::nullopt};
}

void put(std::vector<std::byte>& out, std::uint64_t v, int bytes) {
    for (int b = 0; b < bytes; ++b) out.push_back(static_cast<std::byte>((v >> (8 * b)) & 0xffU));
}

std::uint64_t get(std::span<const std::byte> in, std::size_t at, int bytes) {
    std::uint64_t v = 0;
    for (int b = 0; b < bytes; ++b) {
        v |= static_cast<std::uint64_t>(std::to_integer<unsigned>(in[at + b])) << (8 * b);
    }
    return v;
}

// The container around a payload, written atomically (files.hpp).
std::expected<void, base::Error> write_file(const std::filesystem::path& path, Content kind,
                                            const Writer& payload) {
    const std::vector<std::byte>& body = payload.bytes();
    std::vector<std::byte> out;
    out.reserve(header_size + body.size() + 8);
    for (const char c : magic) out.push_back(static_cast<std::byte>(c));
    put(out, format_version, 4);
    put(out, static_cast<std::uint32_t>(kind), 4);
    put(out, body.size(), 8);
    out.insert(out.end(), body.begin(), body.end());
    put(out, detail::fnv1a(body), 8);

    std::filesystem::path temporary = path;
    temporary += ".tmp";
    {
        std::ofstream f(temporary, std::ios::binary | std::ios::trunc);
        if (!f) return std::unexpected(unwritable(path, "cannot create the temporary file"));
        f.write(reinterpret_cast<const char*>(out.data()), static_cast<std::streamsize>(out.size()));
        f.flush();
        if (!f) {
            f.close();
            std::error_code ignored;
            std::filesystem::remove(temporary, ignored);
            return std::unexpected(unwritable(path, "writing the temporary file failed"));
        }
    }
    std::error_code ec;
    std::filesystem::rename(temporary, path, ec);
    if (ec) {
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        return std::unexpected(unwritable(path, "rename failed: " + ec.message()));
    }
    return {};
}

// The payload of a container of this kind, checked (files.hpp).
std::expected<std::vector<std::byte>, base::Error> read_file(const std::filesystem::path& path,
                                                             Content kind) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return std::unexpected(invalid("cannot open '" + path.string() + "'"));
    std::vector<char> raw((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    std::vector<std::byte> bytes(raw.size());
    for (std::size_t i = 0; i < raw.size(); ++i) bytes[i] = static_cast<std::byte>(raw[i]);
    const std::span<const std::byte> in = bytes;
    if (in.size() < header_size + 8) return std::unexpected(corrupt("shorter than a header"));
    for (std::size_t i = 0; i < magic.size(); ++i) {
        if (std::to_integer<char>(in[i]) != magic[i]) {
            return std::unexpected(invalid("'" + path.string() + "' is not a NiTCAD file"));
        }
    }
    if (const auto v = get(in, 8, 4); v != format_version) {
        return std::unexpected(corrupt("format version " + std::to_string(v) + ", this build reads " +
                                       std::to_string(format_version)));
    }
    if (const auto k = get(in, 12, 4); k != static_cast<std::uint32_t>(kind)) {
        return std::unexpected(corrupt("content kind " + std::to_string(k) + ", expected " +
                                       std::to_string(static_cast<std::uint32_t>(kind))));
    }
    const std::uint64_t length = get(in, 16, 8);
    if (length != in.size() - header_size - 8) return std::unexpected(corrupt("length mismatch"));
    const auto body = in.subspan(header_size, static_cast<std::size_t>(length));
    if (detail::fnv1a(body) != get(in, header_size + body.size(), 8)) {
        return std::unexpected(corrupt("checksum mismatch"));
    }
    return std::vector<std::byte>(body.begin(), body.end());
}

// --- mesh --------------------------------------------------------------------------------------

namespace mesh_key {
constexpr std::uint16_t dimension = 1, points = 2, volumes = 3, endpoints = 4, lengths = 5,
                        couplings = 6, patch = 7;
constexpr std::uint16_t patch_name = 1, patch_nodes = 2, patch_areas = 3;
}  // namespace mesh_key

Writer write_mesh(const mesh::Mesh& m) {
    Writer w;
    w.u64(mesh_key::dimension, static_cast<std::uint64_t>(m.dimension()));
    std::vector<double> xyz;
    xyz.reserve(3 * m.node_count());
    for (const mesh::Point& p : m.points()) xyz.insert(xyz.end(), p.begin(), p.end());
    w.f64s(mesh_key::points, xyz);
    w.f64s(mesh_key::volumes, m.volumes());
    std::vector<std::int64_t> ends;
    std::vector<double> lengths, couplings;
    for (const mesh::Edge& e : m.edges()) {
        ends.push_back(e.first);
        ends.push_back(e.second);
        lengths.push_back(e.length);
        couplings.push_back(e.coupling_area);
    }
    w.i64s(mesh_key::endpoints, ends);
    w.f64s(mesh_key::lengths, lengths);
    w.f64s(mesh_key::couplings, couplings);
    for (const mesh::BoundaryPatch& b : m.boundary()) {
        Writer p;
        p.text(mesh_key::patch_name, b.name);
        std::vector<std::int64_t> nodes(b.nodes.begin(), b.nodes.end());
        p.i64s(mesh_key::patch_nodes, nodes);
        p.f64s(mesh_key::patch_areas, b.areas);
        w.record(mesh_key::patch, p);
    }
    return w;
}

std::expected<mesh::Mesh, base::Error> read_mesh(const Reader& r) {
    NITCAD_IO_TRY(dimension, r.u64(mesh_key::dimension));
    NITCAD_IO_TRY(xyz, r.f64s(mesh_key::points));
    NITCAD_IO_TRY(volumes, r.f64s(mesh_key::volumes));
    NITCAD_IO_TRY(ends, r.i64s(mesh_key::endpoints));
    NITCAD_IO_TRY(lengths, r.f64s(mesh_key::lengths));
    NITCAD_IO_TRY(couplings, r.f64s(mesh_key::couplings));
    NITCAD_IO_TRY(patches, r.records(mesh_key::patch));
    if (xyz.size() % 3 != 0 || ends.size() != 2 * lengths.size() ||
        couplings.size() != lengths.size()) {
        return std::unexpected(corrupt("mesh array sizes"));
    }
    std::vector<mesh::Point> points(xyz.size() / 3);
    for (std::size_t i = 0; i < points.size(); ++i) {
        points[i] = {xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]};
    }
    const auto node_id = [](std::int64_t v) -> std::optional<mesh::NodeId> {
        if (v < 0 || v > std::numeric_limits<mesh::NodeId>::max()) return std::nullopt;
        return static_cast<mesh::NodeId>(v);
    };
    std::vector<mesh::Edge> edges(lengths.size());
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto a = node_id(ends[2 * e]);
        const auto b = node_id(ends[2 * e + 1]);
        if (!a || !b) return std::unexpected(corrupt("edge endpoint out of range"));
        edges[e] = {*a, *b, lengths[e], couplings[e]};
    }
    std::vector<mesh::BoundaryPatch> boundary;
    for (const Reader& p : patches) {
        NITCAD_IO_TRY(name, p.text(mesh_key::patch_name));
        NITCAD_IO_TRY(nodes, p.i64s(mesh_key::patch_nodes));
        NITCAD_IO_TRY(areas, p.f64s(mesh_key::patch_areas));
        mesh::BoundaryPatch b{std::move(name), {}, std::move(areas)};
        for (const std::int64_t v : nodes) {
            const auto id = node_id(v);
            if (!id) return std::unexpected(corrupt("boundary node out of range"));
            b.nodes.push_back(*id);
        }
        boundary.push_back(std::move(b));
    }
    return mesh::Mesh::from_parts(static_cast<int>(dimension), std::move(points),
                                  std::move(volumes), std::move(edges), std::move(boundary));
}

// --- results -----------------------------------------------------------------------------------

namespace fields_key {
constexpr std::uint16_t potential = 1, n = 2, p = 3, temperature = 4;
}
namespace bands_key {
constexpr std::uint16_t conduction = 1, valence = 2, electron = 3, hole = 4;
}
namespace convergence_key {
constexpr std::uint16_t iteration = 1, update = 2, residual = 3, converged = 4, pivot_ratio = 5;
}
namespace paths_key {
constexpr std::uint16_t start = 1, length = 2, field = 3, generation = 4, reached = 5;
}
namespace point_key {
constexpr std::uint16_t bias = 1, fields = 2, bands = 3, current = 4, resolution = 5,
                        gate_charge = 6, trap_charge = 7, edge_n = 8, edge_p = 9,
                        convergence = 10, paths = 11, relocations = 12, thermal_bias = 13,
                        thermal_heat = 14, substeps = 15;
}
namespace run_key {
constexpr std::uint16_t identity = 1, setting = 2, name = 1, value = 2;
}
namespace error_key {
constexpr std::uint16_t code = 1, message = 2, index = 3, value = 4;
}
namespace sweep_key {
constexpr std::uint16_t run = 1, point = 2, stopped = 3, unfinished = 4;
}
namespace checkpoint_key {
constexpr std::uint16_t mesh_digest = 1, run_identity = 2, bias = 3, thermal_bias = 4,
                        fields = 5, nodes = 6;
}

Writer write_fields(const results::NodeFields& f) {
    Writer w;
    w.f64s(fields_key::potential, f.potential_V);
    w.f64s(fields_key::n, f.n_cm3);
    w.f64s(fields_key::p, f.p_cm3);
    w.f64s(fields_key::temperature, f.temperature_K);
    return w;
}

std::expected<results::NodeFields, base::Error> read_fields(const Reader& r) {
    results::NodeFields f;
    NITCAD_IO_TRY(potential, r.f64s_or_empty(fields_key::potential));
    NITCAD_IO_TRY(n, r.f64s_or_empty(fields_key::n));
    NITCAD_IO_TRY(p, r.f64s_or_empty(fields_key::p));
    NITCAD_IO_TRY(temperature, r.f64s_or_empty(fields_key::temperature));
    f.potential_V = std::move(potential);
    f.n_cm3 = std::move(n);
    f.p_cm3 = std::move(p);
    f.temperature_K = std::move(temperature);
    return f;
}

Writer write_convergence(const results::ConvergenceRecord& c) {
    Writer w;
    std::vector<std::int64_t> iteration;
    std::vector<double> update, residual;
    for (const results::IterationRecord& i : c.iterations) {
        iteration.push_back(i.iteration);
        update.push_back(i.update);
        residual.push_back(i.residual);
    }
    w.i64s(convergence_key::iteration, iteration);
    w.f64s(convergence_key::update, update);
    w.f64s(convergence_key::residual, residual);
    w.u64(convergence_key::converged, c.converged ? 1 : 0);
    if (c.smallest_pivot_ratio) w.f64(convergence_key::pivot_ratio, *c.smallest_pivot_ratio);
    return w;
}

std::expected<results::ConvergenceRecord, base::Error> read_convergence(const Reader& r) {
    results::ConvergenceRecord c;
    NITCAD_IO_TRY(iteration, r.i64s(convergence_key::iteration));
    NITCAD_IO_TRY(update, r.f64s(convergence_key::update));
    NITCAD_IO_TRY(residual, r.f64s(convergence_key::residual));
    NITCAD_IO_TRY(converged, r.u64(convergence_key::converged));
    if (update.size() != iteration.size() || residual.size() != iteration.size()) {
        return std::unexpected(corrupt("convergence record sizes"));
    }
    for (std::size_t i = 0; i < iteration.size(); ++i) {
        c.iterations.push_back({static_cast<int>(iteration[i]), update[i], residual[i]});
    }
    c.converged = converged != 0;
    if (r.has(convergence_key::pivot_ratio)) {
        NITCAD_IO_TRY(ratio, r.f64(convergence_key::pivot_ratio));
        c.smallest_pivot_ratio = ratio;
    }
    return c;
}

Writer write_point(const results::BiasPoint& p) {
    Writer w;
    w.f64s(point_key::bias, p.bias_V);
    w.record(point_key::fields, write_fields(p.fields));
    Writer bands;
    bands.f64s(bands_key::conduction, p.bands.conduction_eV);
    bands.f64s(bands_key::valence, p.bands.valence_eV);
    bands.f64s(bands_key::electron, p.bands.electron_fermi_eV);
    bands.f64s(bands_key::hole, p.bands.hole_fermi_eV);
    w.record(point_key::bands, bands);
    w.f64s(point_key::current, p.terminal_current);
    w.f64s(point_key::resolution, p.terminal_current_resolution);
    w.f64s(point_key::gate_charge, p.gate_charge);
    w.f64s(point_key::trap_charge, p.interface_trap_charge);
    w.f64s(point_key::edge_n, p.edge_current_n);
    w.f64s(point_key::edge_p, p.edge_current_p);
    w.record(point_key::convergence, write_convergence(p.convergence));
    Writer paths;
    std::vector<std::int64_t> start, reached;
    std::vector<double> length, field, generation;
    for (const results::TunnelPath& t : p.tunnel_paths) {
        start.push_back(static_cast<std::int64_t>(t.start_node));
        length.push_back(t.length_cm);
        field.push_back(t.field_V_per_cm);
        generation.push_back(t.generation_cm3_s);
        reached.push_back(t.reached ? 1 : 0);
    }
    paths.i64s(paths_key::start, start);
    paths.f64s(paths_key::length, length);
    paths.f64s(paths_key::field, field);
    paths.f64s(paths_key::generation, generation);
    paths.i64s(paths_key::reached, reached);
    w.record(point_key::paths, paths);
    w.i64(point_key::relocations, p.path_relocations);
    w.f64s(point_key::thermal_bias, p.thermal_bias_K);
    w.f64s(point_key::thermal_heat, p.thermal_contact_heat);
    w.i64(point_key::substeps, p.bias_substeps);
    return w;
}

std::expected<results::BiasPoint, base::Error> read_point(const Reader& r) {
    results::BiasPoint p;
    NITCAD_IO_TRY(bias, r.f64s(point_key::bias));
    p.bias_V = std::move(bias);
    NITCAD_IO_TRY(fields_record, r.record(point_key::fields));
    NITCAD_IO_TRY(fields, read_fields(fields_record));
    p.fields = std::move(fields);
    if (r.has(point_key::bands)) {
        NITCAD_IO_TRY(b, r.record(point_key::bands));
        NITCAD_IO_TRY(conduction, b.f64s_or_empty(bands_key::conduction));
        NITCAD_IO_TRY(valence, b.f64s_or_empty(bands_key::valence));
        NITCAD_IO_TRY(electron, b.f64s_or_empty(bands_key::electron));
        NITCAD_IO_TRY(hole, b.f64s_or_empty(bands_key::hole));
        p.bands = {std::move(conduction), std::move(valence), std::move(electron), std::move(hole)};
    }
    NITCAD_IO_TRY(current, r.f64s_or_empty(point_key::current));
    NITCAD_IO_TRY(resolution, r.f64s_or_empty(point_key::resolution));
    NITCAD_IO_TRY(gate_charge, r.f64s_or_empty(point_key::gate_charge));
    NITCAD_IO_TRY(trap_charge, r.f64s_or_empty(point_key::trap_charge));
    NITCAD_IO_TRY(edge_n, r.f64s_or_empty(point_key::edge_n));
    NITCAD_IO_TRY(edge_p, r.f64s_or_empty(point_key::edge_p));
    p.terminal_current = std::move(current);
    p.terminal_current_resolution = std::move(resolution);
    p.gate_charge = std::move(gate_charge);
    p.interface_trap_charge = std::move(trap_charge);
    p.edge_current_n = std::move(edge_n);
    p.edge_current_p = std::move(edge_p);
    if (r.has(point_key::convergence)) {
        NITCAD_IO_TRY(c, r.record(point_key::convergence));
        NITCAD_IO_TRY(convergence, read_convergence(c));
        p.convergence = std::move(convergence);
    }
    if (r.has(point_key::paths)) {
        NITCAD_IO_TRY(t, r.record(point_key::paths));
        NITCAD_IO_TRY(start, t.i64s(paths_key::start));
        NITCAD_IO_TRY(length, t.f64s(paths_key::length));
        NITCAD_IO_TRY(field, t.f64s(paths_key::field));
        NITCAD_IO_TRY(generation, t.f64s(paths_key::generation));
        NITCAD_IO_TRY(reached, t.i64s(paths_key::reached));
        if (length.size() != start.size() || field.size() != start.size() ||
            generation.size() != start.size() || reached.size() != start.size()) {
            return std::unexpected(corrupt("tunnel path sizes"));
        }
        for (std::size_t i = 0; i < start.size(); ++i) {
            if (start[i] < 0) return std::unexpected(corrupt("tunnel path start"));
            p.tunnel_paths.push_back({static_cast<std::size_t>(start[i]), length[i], field[i],
                                      generation[i], reached[i] != 0});
        }
    }
    NITCAD_IO_TRY(relocations, r.i64_or(point_key::relocations, 0));
    NITCAD_IO_TRY(thermal_bias, r.f64s_or_empty(point_key::thermal_bias));
    NITCAD_IO_TRY(thermal_heat, r.f64s_or_empty(point_key::thermal_heat));
    NITCAD_IO_TRY(substeps, r.i64_or(point_key::substeps, 0));
    p.path_relocations = static_cast<int>(relocations);
    p.thermal_bias_K = std::move(thermal_bias);
    p.thermal_contact_heat = std::move(thermal_heat);
    p.bias_substeps = static_cast<int>(substeps);
    return p;
}

Writer write_run(const results::RunRecord& run) {
    Writer w;
    w.u64(run_key::identity, run.input_identity);
    for (const auto& [name, value] : run.settings) {
        Writer s;
        s.text(run_key::name, name);
        s.f64(run_key::value, value);
        w.record(run_key::setting, s);
    }
    return w;
}

std::expected<results::RunRecord, base::Error> read_run(const Reader& r) {
    results::RunRecord run;
    NITCAD_IO_TRY(identity, r.u64(run_key::identity));
    run.input_identity = identity;
    NITCAD_IO_TRY(settings, r.records(run_key::setting));
    for (const Reader& s : settings) {
        NITCAD_IO_TRY(name, s.text(run_key::name));
        NITCAD_IO_TRY(value, s.f64(run_key::value));
        run.settings.emplace_back(std::move(name), value);
    }
    return run;
}

Writer write_error(const base::Error& e) {
    Writer w;
    w.u64(error_key::code, static_cast<std::uint64_t>(e.code));
    w.text(error_key::message, e.message);
    if (e.context && e.context->index) w.u64(error_key::index, *e.context->index);
    if (e.context && e.context->value) w.f64(error_key::value, *e.context->value);
    return w;
}

std::expected<base::Error, base::Error> read_error(const Reader& r) {
    NITCAD_IO_TRY(code, r.u64(error_key::code));
    if (code > static_cast<std::uint64_t>(base::ErrorCode::resource_exhausted)) {
        return std::unexpected(corrupt("error code"));
    }
    NITCAD_IO_TRY(message, r.text(error_key::message));
    base::Error e{static_cast<base::ErrorCode>(code), std::move(message), std::nullopt};
    if (r.has(error_key::index) || r.has(error_key::value)) {
        base::ErrorContext c;
        if (r.has(error_key::index)) {
            NITCAD_IO_TRY(index, r.u64(error_key::index));
            c.index = static_cast<std::size_t>(index);
        }
        if (r.has(error_key::value)) {
            NITCAD_IO_TRY(value, r.f64(error_key::value));
            c.value = value;
        }
        e.context = c;
    }
    return e;
}

void write_digest_reals(std::uint64_t& h, std::span<const double> v) {
    for (const double d : v) {
        std::array<std::byte, 8> b{};
        const auto bits = std::bit_cast<std::uint64_t>(d);
        for (int i = 0; i < 8; ++i) b[i] = static_cast<std::byte>((bits >> (8 * i)) & 0xffU);
        h = detail::fnv1a(b, h);
    }
}

void write_digest_integer(std::uint64_t& h, std::uint64_t v) {
    std::array<std::byte, 8> b{};
    for (int i = 0; i < 8; ++i) b[i] = static_cast<std::byte>((v >> (8 * i)) & 0xffU);
    h = detail::fnv1a(b, h);
}

}  // namespace

std::uint64_t mesh_digest(const mesh::Mesh& m) {
    std::uint64_t h = 0xcbf29ce484222325ULL;
    write_digest_integer(h, static_cast<std::uint64_t>(m.dimension()));
    write_digest_integer(h, m.node_count());
    for (const mesh::Point& p : m.points()) write_digest_reals(h, p);
    write_digest_reals(h, m.volumes());
    write_digest_integer(h, m.edges().size());
    for (const mesh::Edge& e : m.edges()) {
        write_digest_integer(h, static_cast<std::uint64_t>(e.first));
        write_digest_integer(h, static_cast<std::uint64_t>(e.second));
        write_digest_reals(h, std::array{e.length, e.coupling_area});
    }
    write_digest_integer(h, m.boundary().size());
    for (const mesh::BoundaryPatch& b : m.boundary()) {
        h = detail::fnv1a(std::as_bytes(std::span(b.name.data(), b.name.size())), h);
        write_digest_integer(h, b.nodes.size());
        for (const mesh::NodeId n : b.nodes) write_digest_integer(h, static_cast<std::uint64_t>(n));
        write_digest_reals(h, b.areas);
    }
    return h;
}

std::expected<void, base::Error> save_mesh(const std::filesystem::path& path,
                                           const mesh::Mesh& m) {
    return write_file(path, Content::mesh, write_mesh(m));
}

std::expected<mesh::Mesh, base::Error> load_mesh(const std::filesystem::path& path) {
    NITCAD_IO_TRY(bytes, read_file(path, Content::mesh));
    NITCAD_IO_TRY(r, Reader::parse(bytes));
    return read_mesh(r);
}

std::expected<void, base::Error> save_sweep(const std::filesystem::path& path,
                                            const results::Sweep& sweep) {
    Writer w;
    w.record(sweep_key::run, write_run(sweep.run));
    for (const results::BiasPoint& p : sweep.points) w.record(sweep_key::point, write_point(p));
    if (sweep.stopped) w.record(sweep_key::stopped, write_error(*sweep.stopped));
    if (sweep.unfinished) w.record(sweep_key::unfinished, write_convergence(*sweep.unfinished));
    return write_file(path, Content::sweep, w);
}

std::expected<results::Sweep, base::Error> load_sweep(const std::filesystem::path& path) {
    NITCAD_IO_TRY(bytes, read_file(path, Content::sweep));
    NITCAD_IO_TRY(r, Reader::parse(bytes));
    results::Sweep sweep;
    NITCAD_IO_TRY(run_record, r.record(sweep_key::run));
    NITCAD_IO_TRY(run, read_run(run_record));
    sweep.run = std::move(run);
    NITCAD_IO_TRY(points, r.records(sweep_key::point));
    for (const Reader& p : points) {
        NITCAD_IO_TRY(point, read_point(p));
        sweep.points.push_back(std::move(point));
    }
    if (r.has(sweep_key::stopped)) {
        NITCAD_IO_TRY(e, r.record(sweep_key::stopped));
        NITCAD_IO_TRY(stopped, read_error(e));
        sweep.stopped = std::move(stopped);
    }
    if (r.has(sweep_key::unfinished)) {
        NITCAD_IO_TRY(c, r.record(sweep_key::unfinished));
        NITCAD_IO_TRY(unfinished, read_convergence(c));
        sweep.unfinished = std::move(unfinished);
    }
    return sweep;
}

std::expected<void, base::Error> save_checkpoint(const std::filesystem::path& path,
                                                 const mesh::Mesh& m,
                                                 const results::RunRecord& run,
                                                 const results::BiasPoint& point) {
    Writer w;
    w.u64(checkpoint_key::mesh_digest, mesh_digest(m));
    w.u64(checkpoint_key::nodes, m.node_count());
    w.u64(checkpoint_key::run_identity, run.input_identity);
    w.f64s(checkpoint_key::bias, point.bias_V);
    w.f64s(checkpoint_key::thermal_bias, point.thermal_bias_K);
    w.record(checkpoint_key::fields, write_fields(point.fields));
    return write_file(path, Content::checkpoint, w);
}

std::expected<Checkpoint, base::Error> load_checkpoint(const std::filesystem::path& path,
                                                       const mesh::Mesh& m) {
    NITCAD_IO_TRY(bytes, read_file(path, Content::checkpoint));
    NITCAD_IO_TRY(r, Reader::parse(bytes));
    Checkpoint c;
    NITCAD_IO_TRY(digest, r.u64(checkpoint_key::mesh_digest));
    if (digest != mesh_digest(m)) {
        return std::unexpected(invalid("checkpoint '" + path.string() +
                                       "' was written for another mesh"));
    }
    c.mesh_digest = digest;
    NITCAD_IO_TRY(identity, r.u64(checkpoint_key::run_identity));
    c.run_identity = identity;
    NITCAD_IO_TRY(bias, r.f64s(checkpoint_key::bias));
    c.bias_V = std::move(bias);
    NITCAD_IO_TRY(thermal, r.f64s_or_empty(checkpoint_key::thermal_bias));
    c.thermal_bias_K = std::move(thermal);
    NITCAD_IO_TRY(fields_record, r.record(checkpoint_key::fields));
    NITCAD_IO_TRY(fields, read_fields(fields_record));
    const std::size_t n = m.node_count();
    const auto sized = [n](const std::vector<double>& v) { return v.size() == n; };
    if (!sized(fields.potential_V) || !sized(fields.n_cm3) || !sized(fields.p_cm3) ||
        !(fields.temperature_K.empty() || sized(fields.temperature_K))) {
        return std::unexpected(corrupt("checkpoint fields do not have one value per node"));
    }
    c.fields = std::move(fields);
    return c;
}

}  // namespace NiTCAD::io

#undef NITCAD_IO_TRY
