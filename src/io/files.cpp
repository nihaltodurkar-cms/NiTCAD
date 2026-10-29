#include "NiTCAD/io/files.hpp"

#include <array>
#include <bit>
#include <complex>
#include <cstddef>
#include <fstream>
#include <iterator>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <system_error>
#include <utility>

#include "container.hpp"
#include "record.hpp"

namespace NiTCAD::io {

namespace detail {

base::Error invalid(std::string message) {
    return {base::ErrorCode::invalid_input, std::move(message), std::nullopt};
}

}  // namespace detail

namespace {

using detail::corrupt;
using detail::invalid;
using detail::Reader;
using detail::Writer;

constexpr std::array<char, 8> magic{'N', 'i', 'T', 'C', 'A', 'D', 'i', 'o'};
constexpr std::size_t header_size = 8 + 4 + 4 + 8;

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

}  // namespace

namespace detail {

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
        return std::unexpected(corrupt("format version " + std::to_string(v) +
                                       ", this build reads " +
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
                        couplings = 6, patch = 7, origin = 8,  // origin: absent in older files
                        // A generated mesh (M14), absent in plain mesh files:
                        cells = 9, facets = 10, facet_patch = 11, node_region = 12, report = 13,
                        certificate = 14;
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
    w.u64(mesh_key::origin, static_cast<std::uint64_t>(m.origin()));
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
    NITCAD_IO_TRY(origin, r.u64_or(mesh_key::origin, 0));  // a file without it: a general graph
    if (origin > static_cast<std::uint64_t>(mesh::MeshOrigin::tensor_grid)) {
        return std::unexpected(corrupt("unknown mesh origin"));
    }
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
                                  std::move(volumes), std::move(edges), std::move(boundary),
                                  static_cast<mesh::MeshOrigin>(origin));
}

namespace certificate_key {
constexpr std::uint16_t digests = 1, counts = 2, deviations = 3;
}  // namespace certificate_key

Writer write_certificate(const mesh::InterfaceCertificate& c) {
    Writer cw;
    cw.i64s(certificate_key::digests,
            std::vector<std::int64_t>{std::bit_cast<std::int64_t>(c.mesh_digest),
                                      std::bit_cast<std::int64_t>(c.region_digest),
                                      std::bit_cast<std::int64_t>(c.interfaces_digest)});
    cw.i64s(certificate_key::counts,
            std::vector<std::int64_t>{static_cast<std::int64_t>(c.report.cross_faces),
                                      static_cast<std::int64_t>(c.report.corner_faces)});
    cw.f64s(certificate_key::deviations,
            std::vector<double>{c.report.max_flat_deviation, c.report.max_corner_deviation});
    return cw;
}

std::expected<mesh::InterfaceCertificate, base::Error> read_certificate(const Reader& cr) {
    NITCAD_IO_TRY(digests, cr.i64s(certificate_key::digests));
    NITCAD_IO_TRY(ccounts, cr.i64s(certificate_key::counts));
    NITCAD_IO_TRY(deviations, cr.f64s(certificate_key::deviations));
    if (digests.size() != 3 || ccounts.size() != 2 || deviations.size() != 2) {
        return std::unexpected(corrupt("certificate sizes"));
    }
    return mesh::InterfaceCertificate{.mesh_digest = std::bit_cast<std::uint64_t>(digests[0]),
                                      .region_digest = std::bit_cast<std::uint64_t>(digests[1]),
                                      .interfaces_digest = std::bit_cast<std::uint64_t>(digests[2]),
                                      .report = {static_cast<std::size_t>(ccounts[0]),
                                                 static_cast<std::size_t>(ccounts[1]), deviations[0],
                                                 deviations[1]}};
}

}  // namespace detail

namespace {

using detail::read_file;
using detail::read_mesh;
using detail::write_file;
using detail::write_mesh;
namespace mesh_key = detail::mesh_key;

// --- results -----------------------------------------------------------------------------------

namespace fields_key {
constexpr std::uint16_t potential = 1, n = 2, p = 3, temperature = 4;
}
namespace rates_key {  // Unit 45
constexpr std::uint16_t srh = 1, auger = 2, radiative = 3, impact = 4, interface_n = 5, interface_p = 6,
                        band_to_band_n = 7, band_to_band_p = 8, net_n = 9, net_p = 10;
}
namespace bands_key {
constexpr std::uint16_t conduction = 1, valence = 2, electron = 3, hole = 4;
}
namespace convergence_key {
constexpr std::uint16_t iteration = 1, update = 2, residual = 3, converged = 4, pivot_ratio = 5,
                        step = 6;
}
namespace paths_key {
constexpr std::uint16_t start = 1, length = 2, field = 3, generation = 4, reached = 5;
}
namespace point_key {
constexpr std::uint16_t bias = 1, fields = 2, bands = 3, current = 4, resolution = 5,
                        gate_charge = 6, trap_charge = 7, edge_n = 8, edge_p = 9,
                        convergence = 10, paths = 11, relocations = 12, thermal_bias = 13,
                        thermal_heat = 14, substeps = 15, weighted = 16,
                        weighted_resolution = 17;  // Unit 32: absent when off
constexpr std::uint16_t rates = 18;                // Unit 45: absent when not stored
constexpr std::uint16_t tunnelling_ecb = 19, tunnelling_hvb = 20,
                        tunnelling_evb = 21;       // Unit 51: absent without gate tunnelling
}
namespace run_key {
constexpr std::uint16_t identity = 1, setting = 2, name = 1, value = 2;
constexpr std::uint16_t provenance = 3;   // Unit 42; absent in older files
constexpr std::uint16_t mesh_digest = 4;  // Unit 45; absent in older files (0: unknown)
}
namespace provenance_key {
constexpr std::uint16_t version = 1, revision = 2, dirty = 3, build_type = 4, compiler = 5,
                        eigen_version = 6, eigen_cache = 7, mkl_version = 8, os = 9, cpu = 10,
                        logical_cpus = 11, started = 12, wall = 13;
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
constexpr std::uint16_t provenance = 7;  // Unit 42; absent in older files
}
// Transient runs (Unit 42).
namespace transient_key {
constexpr std::uint16_t run = 1, point = 2, snapshot = 3, rejected = 4, retraced = 5, lagged = 6,
                        stopped = 7, unfinished = 8;
}
namespace time_key {
constexpr std::uint16_t time = 1, step = 2, order = 3, error_ratio = 4, bias = 5, current = 6,
                        conduction = 7, displacement = 8, contact_charge = 9, trap_charge = 10,
                        thermal_heat = 11, convergence = 12;
}
namespace snapshot_key {
constexpr std::uint16_t time = 1, point = 2, fields = 3, bands = 4, edge_n = 5, edge_p = 6,
                        occupancy = 7;
constexpr std::uint16_t rates = 8;  // Unit 45: absent when not stored
}
// Small-signal runs (Unit 42). Complex values are stored as (real, imaginary) pairs of doubles.
namespace small_signal_key {
constexpr std::uint16_t run = 1, contacts = 2, frequencies = 3, point = 4, stopped = 5,
                        unfinished = 6;
}
namespace ac_point_key {
// admittance and resolution repeat once per frequency, in frequency order.
constexpr std::uint16_t dc = 1, admittance = 2, pivot_ratio = 3, backward_error = 4,
                        resolution = 5, fields = 6;
}
namespace ac_fields_key {
constexpr std::uint16_t frequency = 1, contact = 2, potential = 3, n = 4, p = 5, temperature = 6;
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

// The rates of a point or a snapshot (Unit 45): every vector, empty ones as absent.
Writer write_rates(const results::RateFields& r) {
    namespace k = rates_key;
    Writer w;
    w.f64s(k::srh, r.srh);
    w.f64s(k::auger, r.auger);
    w.f64s(k::radiative, r.radiative);
    w.f64s(k::impact, r.impact);
    w.f64s(k::interface_n, r.interface_n);
    w.f64s(k::interface_p, r.interface_p);
    w.f64s(k::band_to_band_n, r.band_to_band_n);
    w.f64s(k::band_to_band_p, r.band_to_band_p);
    w.f64s(k::net_n, r.net_n);
    w.f64s(k::net_p, r.net_p);
    return w;
}

std::expected<results::RateFields, base::Error> read_rates(const Reader& r) {
    namespace k = rates_key;
    results::RateFields out;
    const std::array<std::pair<std::uint16_t, std::vector<double>*>, 10> slots{
        {{k::srh, &out.srh}, {k::auger, &out.auger}, {k::radiative, &out.radiative},
         {k::impact, &out.impact}, {k::interface_n, &out.interface_n}, {k::interface_p, &out.interface_p},
         {k::band_to_band_n, &out.band_to_band_n}, {k::band_to_band_p, &out.band_to_band_p},
         {k::net_n, &out.net_n}, {k::net_p, &out.net_p}}};
    for (const auto& [key, target] : slots) {
        NITCAD_IO_TRY(values, r.f64s_or_empty(key));
        *target = std::move(values);
    }
    return out;
}

Writer write_convergence(const results::ConvergenceRecord& c) {
    Writer w;
    std::vector<std::int64_t> iteration;
    std::vector<double> update, residual, step;
    for (const results::IterationRecord& i : c.iterations) {
        iteration.push_back(i.iteration);
        update.push_back(i.update);
        residual.push_back(i.residual);
        step.push_back(i.step);
    }
    w.i64s(convergence_key::iteration, iteration);
    w.f64s(convergence_key::update, update);
    w.f64s(convergence_key::residual, residual);
    w.f64s(convergence_key::step, step);
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
    // The fraction of each step taken (Unit 30); absent in files written before it: full steps.
    NITCAD_IO_TRY(step, r.f64s_or_empty(convergence_key::step));
    if (step.empty()) step.assign(iteration.size(), 1.0);
    if (update.size() != iteration.size() || residual.size() != iteration.size() ||
        step.size() != iteration.size()) {
        return std::unexpected(corrupt("convergence record sizes"));
    }
    for (std::size_t i = 0; i < iteration.size(); ++i) {
        c.iterations.push_back({static_cast<int>(iteration[i]), update[i], residual[i], step[i]});
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
    if (!p.weighted_current.empty()) {
        w.f64s(point_key::weighted, p.weighted_current);
        w.f64s(point_key::weighted_resolution, p.weighted_current_resolution);
    }
    if (!p.rates.net_n.empty()) w.record(point_key::rates, write_rates(p.rates));
    if (!p.tunnelling_ecb.empty()) {
        w.f64s(point_key::tunnelling_ecb, p.tunnelling_ecb);
        w.f64s(point_key::tunnelling_hvb, p.tunnelling_hvb);
        w.f64s(point_key::tunnelling_evb, p.tunnelling_evb);
    }
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
    NITCAD_IO_TRY(weighted, r.f64s_or_empty(point_key::weighted));
    NITCAD_IO_TRY(weighted_resolution, r.f64s_or_empty(point_key::weighted_resolution));
    p.weighted_current = std::move(weighted);
    p.weighted_current_resolution = std::move(weighted_resolution);
    if (r.has(point_key::rates)) {
        NITCAD_IO_TRY(rr, r.record(point_key::rates));
        NITCAD_IO_TRY(rates, read_rates(rr));
        p.rates = std::move(rates);
    }
    NITCAD_IO_TRY(ecb, r.f64s_or_empty(point_key::tunnelling_ecb));
    NITCAD_IO_TRY(hvb, r.f64s_or_empty(point_key::tunnelling_hvb));
    NITCAD_IO_TRY(evb, r.f64s_or_empty(point_key::tunnelling_evb));
    p.tunnelling_ecb = std::move(ecb);
    p.tunnelling_hvb = std::move(hvb);
    p.tunnelling_evb = std::move(evb);
    return p;
}

Writer write_provenance(const results::Provenance& p) {
    namespace k = provenance_key;
    Writer w;
    w.text(k::version, p.version);
    w.text(k::revision, p.revision);
    w.u64(k::dirty, p.dirty ? 1 : 0);
    w.text(k::build_type, p.build_type);
    w.text(k::compiler, p.compiler);
    w.text(k::eigen_version, p.eigen_version);
    w.i64s(k::eigen_cache, p.eigen_cache_bytes);
    w.text(k::mkl_version, p.mkl_version);
    w.text(k::os, p.os);
    w.text(k::cpu, p.cpu);
    w.u64(k::logical_cpus, p.logical_cpus);
    w.text(k::started, p.started_utc);
    w.f64(k::wall, p.wall_s);
    return w;
}

// Every field is optional (IO2): a provenance written by a later revision may lack one.
std::expected<results::Provenance, base::Error> read_provenance(const Reader& r) {
    namespace k = provenance_key;
    results::Provenance p;
    const auto text = [&r](std::uint16_t key, std::string& out) -> std::expected<void, base::Error> {
        if (!r.has(key)) return {};
        NITCAD_IO_TRY(v, r.text(key));
        out = std::move(v);
        return {};
    };
    for (const auto& [key, out] :
         {std::pair{k::version, &p.version}, std::pair{k::revision, &p.revision},
          std::pair{k::build_type, &p.build_type}, std::pair{k::compiler, &p.compiler},
          std::pair{k::eigen_version, &p.eigen_version}, std::pair{k::mkl_version, &p.mkl_version},
          std::pair{k::os, &p.os}, std::pair{k::cpu, &p.cpu},
          std::pair{k::started, &p.started_utc}}) {
        if (auto t = text(key, *out); !t) return std::unexpected(std::move(t.error()));
    }
    NITCAD_IO_TRY(dirty, r.u64_or(k::dirty, 0));
    if (dirty > 1) return std::unexpected(corrupt("provenance dirty flag"));
    p.dirty = dirty == 1;
    if (r.has(k::eigen_cache)) {
        NITCAD_IO_TRY(cache, r.i64s(k::eigen_cache));
        p.eigen_cache_bytes = std::move(cache);
    }
    NITCAD_IO_TRY(cpus, r.u64_or(k::logical_cpus, 0));
    if (cpus > std::numeric_limits<std::uint32_t>::max()) {
        return std::unexpected(corrupt("provenance logical CPU count"));
    }
    p.logical_cpus = static_cast<std::uint32_t>(cpus);
    if (r.has(k::wall)) {
        NITCAD_IO_TRY(wall, r.f64(k::wall));
        p.wall_s = wall;
    }
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
    w.record(run_key::provenance, write_provenance(run.provenance));
    if (run.mesh_digest != 0) w.u64(run_key::mesh_digest, run.mesh_digest);
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
    if (r.has(run_key::provenance)) {
        NITCAD_IO_TRY(p, r.record(run_key::provenance));
        NITCAD_IO_TRY(provenance, read_provenance(p));
        run.provenance = std::move(provenance);
    }
    if (r.has(run_key::mesh_digest)) {
        NITCAD_IO_TRY(digest, r.u64(run_key::mesh_digest));
        run.mesh_digest = digest;
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

// --- transient and small-signal runs (Unit 42) -------------------------------------------------

Writer write_bands(const results::BandDiagram& b) {
    Writer w;
    w.f64s(bands_key::conduction, b.conduction_eV);
    w.f64s(bands_key::valence, b.valence_eV);
    w.f64s(bands_key::electron, b.electron_fermi_eV);
    w.f64s(bands_key::hole, b.hole_fermi_eV);
    return w;
}

std::expected<results::BandDiagram, base::Error> read_bands(const Reader& r) {
    NITCAD_IO_TRY(conduction, r.f64s_or_empty(bands_key::conduction));
    NITCAD_IO_TRY(valence, r.f64s_or_empty(bands_key::valence));
    NITCAD_IO_TRY(electron, r.f64s_or_empty(bands_key::electron));
    NITCAD_IO_TRY(hole, r.f64s_or_empty(bands_key::hole));
    return results::BandDiagram{std::move(conduction), std::move(valence), std::move(electron),
                                std::move(hole)};
}

// A size_t count or index stored as u64; corrupt if it does not fit.
std::expected<std::size_t, base::Error> size_of(std::uint64_t v, const char* what) {
    if (v > std::numeric_limits<std::size_t>::max()) return std::unexpected(corrupt(what));
    return static_cast<std::size_t>(v);
}

std::vector<double> interleaved(std::span<const std::complex<double>> v) {
    std::vector<double> out;
    out.reserve(2 * v.size());
    for (const std::complex<double> z : v) {
        out.push_back(z.real());
        out.push_back(z.imag());
    }
    return out;
}

std::expected<std::vector<std::complex<double>>, base::Error> complexes(
    const std::vector<double>& v) {
    if (v.size() % 2 != 0) return std::unexpected(corrupt("odd length of a complex array"));
    std::vector<std::complex<double>> out(v.size() / 2);
    for (std::size_t i = 0; i < out.size(); ++i) out[i] = {v[2 * i], v[2 * i + 1]};
    return out;
}

std::expected<std::vector<std::complex<double>>, base::Error> read_complexes(const Reader& r,
                                                                             std::uint16_t key) {
    NITCAD_IO_TRY(v, r.f64s_or_empty(key));
    return complexes(v);
}

Writer write_time_point(const results::TimePoint& p) {
    namespace k = time_key;
    Writer w;
    w.f64(k::time, p.time_s);
    w.f64(k::step, p.step_s);
    w.i64(k::order, p.order);
    w.f64(k::error_ratio, p.error_ratio);
    w.f64s(k::bias, p.bias_V);
    w.f64s(k::current, p.terminal_current);
    w.f64s(k::conduction, p.conduction_current);
    w.f64s(k::displacement, p.displacement_current);
    w.f64s(k::contact_charge, p.contact_charge);
    w.f64s(k::trap_charge, p.interface_trap_charge);
    w.f64s(k::thermal_heat, p.thermal_contact_heat);
    w.record(k::convergence, write_convergence(p.convergence));
    return w;
}

std::expected<results::TimePoint, base::Error> read_time_point(const Reader& r) {
    namespace k = time_key;
    NITCAD_IO_TRY(time, r.f64(k::time));
    NITCAD_IO_TRY(step, r.f64(k::step));
    NITCAD_IO_TRY(order, r.i64(k::order));
    if (order < 0 || order > 2) return std::unexpected(corrupt("time step order"));
    NITCAD_IO_TRY(error_ratio, r.f64(k::error_ratio));
    NITCAD_IO_TRY(bias, r.f64s_or_empty(k::bias));
    NITCAD_IO_TRY(current, r.f64s_or_empty(k::current));
    NITCAD_IO_TRY(conduction, r.f64s_or_empty(k::conduction));
    NITCAD_IO_TRY(displacement, r.f64s_or_empty(k::displacement));
    NITCAD_IO_TRY(contact_charge, r.f64s_or_empty(k::contact_charge));
    NITCAD_IO_TRY(trap_charge, r.f64s_or_empty(k::trap_charge));
    NITCAD_IO_TRY(thermal_heat, r.f64s_or_empty(k::thermal_heat));
    NITCAD_IO_TRY(c, r.record(k::convergence));
    NITCAD_IO_TRY(convergence, read_convergence(c));
    return results::TimePoint{.time_s = time,
                              .step_s = step,
                              .order = static_cast<int>(order),
                              .error_ratio = error_ratio,
                              .bias_V = std::move(bias),
                              .terminal_current = std::move(current),
                              .conduction_current = std::move(conduction),
                              .displacement_current = std::move(displacement),
                              .contact_charge = std::move(contact_charge),
                              .interface_trap_charge = std::move(trap_charge),
                              .thermal_contact_heat = std::move(thermal_heat),
                              .convergence = std::move(convergence)};
}

Writer write_snapshot(const results::TransientSnapshot& s) {
    namespace k = snapshot_key;
    Writer w;
    w.f64(k::time, s.time_s);
    w.u64(k::point, s.point);
    w.record(k::fields, write_fields(s.fields));
    w.record(k::bands, write_bands(s.bands));
    w.f64s(k::edge_n, s.edge_current_n);
    w.f64s(k::edge_p, s.edge_current_p);
    w.f64s(k::occupancy, s.trap_occupancy);
    if (!s.rates.net_n.empty()) w.record(k::rates, write_rates(s.rates));
    return w;
}

std::expected<results::TransientSnapshot, base::Error> read_snapshot(const Reader& r,
                                                                     std::size_t points) {
    namespace k = snapshot_key;
    results::TransientSnapshot s;
    NITCAD_IO_TRY(time, r.f64(k::time));
    s.time_s = time;
    NITCAD_IO_TRY(point, r.u64(k::point));
    if (point >= points) return std::unexpected(corrupt("snapshot of a time point that is not there"));
    s.point = static_cast<std::size_t>(point);
    NITCAD_IO_TRY(f, r.record(k::fields));
    NITCAD_IO_TRY(fields, read_fields(f));
    s.fields = std::move(fields);
    NITCAD_IO_TRY(b, r.record(k::bands));
    NITCAD_IO_TRY(bands, read_bands(b));
    s.bands = std::move(bands);
    NITCAD_IO_TRY(edge_n, r.f64s_or_empty(k::edge_n));
    NITCAD_IO_TRY(edge_p, r.f64s_or_empty(k::edge_p));
    NITCAD_IO_TRY(occupancy, r.f64s_or_empty(k::occupancy));
    s.edge_current_n = std::move(edge_n);
    s.edge_current_p = std::move(edge_p);
    s.trap_occupancy = std::move(occupancy);
    if (r.has(k::rates)) {
        NITCAD_IO_TRY(rr, r.record(k::rates));
        NITCAD_IO_TRY(rates, read_rates(rr));
        s.rates = std::move(rates);
    }
    return s;
}

Writer write_ac_fields(const results::SmallSignalFields& f) {
    namespace k = ac_fields_key;
    Writer w;
    w.f64(k::frequency, f.frequency_Hz);
    w.u64(k::contact, f.contact);
    w.f64s(k::potential, interleaved(f.potential));
    w.f64s(k::n, interleaved(f.n_cm3));
    w.f64s(k::p, interleaved(f.p_cm3));
    w.f64s(k::temperature, interleaved(f.temperature_K));
    return w;
}

std::expected<results::SmallSignalFields, base::Error> read_ac_fields(const Reader& r) {
    namespace k = ac_fields_key;
    results::SmallSignalFields f;
    NITCAD_IO_TRY(frequency, r.f64(k::frequency));
    f.frequency_Hz = frequency;
    NITCAD_IO_TRY(contact, r.u64(k::contact));
    NITCAD_IO_TRY(c, size_of(contact, "small-signal field contact"));
    f.contact = c;
    NITCAD_IO_TRY(potential, read_complexes(r, k::potential));
    NITCAD_IO_TRY(n, read_complexes(r, k::n));
    NITCAD_IO_TRY(p, read_complexes(r, k::p));
    NITCAD_IO_TRY(temperature, read_complexes(r, k::temperature));
    f.potential = std::move(potential);
    f.n_cm3 = std::move(n);
    f.p_cm3 = std::move(p);
    f.temperature_K = std::move(temperature);
    return f;
}

// The admittance matrices of all frequencies as one array, in frequency order, each row-major as
// in memory; the resolutions likewise (frequency, then driven contact).
Writer write_ac_point(const results::SmallSignalPoint& p) {
    namespace k = ac_point_key;
    Writer w;
    w.record(k::dc, write_point(p.dc));
    std::vector<double> admittance, resolution;
    for (const auto& y : p.admittance) {
        const std::vector<double> v = interleaved(y);
        admittance.insert(admittance.end(), v.begin(), v.end());
    }
    for (const auto& r : p.resolution) resolution.insert(resolution.end(), r.begin(), r.end());
    w.f64s(k::admittance, admittance);
    w.f64s(k::pivot_ratio, p.pivot_ratio);
    w.f64s(k::backward_error, p.backward_error);
    w.f64s(k::resolution, resolution);
    for (const auto& f : p.fields) w.record(k::fields, write_ac_fields(f));
    return w;
}

// Every point has one contacts x contacts admittance matrix per frequency, so that
// SmallSignal::admittance stays in bounds; the resolution has one value per frequency and driven
// contact, or none.
std::expected<results::SmallSignalPoint, base::Error> read_ac_point(const Reader& r,
                                                                    std::size_t frequencies,
                                                                    std::size_t contacts) {
    namespace k = ac_point_key;
    results::SmallSignalPoint p;
    NITCAD_IO_TRY(dc, r.record(k::dc));
    NITCAD_IO_TRY(point, read_point(dc));
    p.dc = std::move(point);
    NITCAD_IO_TRY(flat, r.f64s_or_empty(k::admittance));
    NITCAD_IO_TRY(admittance, complexes(flat));
    // load_small_signal bounds contacts and frequencies, so these products do not overflow.
    const std::size_t matrix = contacts * contacts;
    if (admittance.size() != frequencies * matrix) {
        return std::unexpected(corrupt("admittances do not fit the frequencies and contacts"));
    }
    for (std::size_t f = 0; f < frequencies; ++f) {
        const auto first = admittance.begin() + static_cast<std::ptrdiff_t>(f * matrix);
        p.admittance.emplace_back(first, first + static_cast<std::ptrdiff_t>(matrix));
    }
    NITCAD_IO_TRY(pivot_ratio, r.f64s_or_empty(k::pivot_ratio));
    NITCAD_IO_TRY(backward_error, r.f64s_or_empty(k::backward_error));
    p.pivot_ratio = std::move(pivot_ratio);
    p.backward_error = std::move(backward_error);
    NITCAD_IO_TRY(resolution, r.f64s_or_empty(k::resolution));
    if (!resolution.empty()) {
        if (resolution.size() != frequencies * contacts) {
            return std::unexpected(corrupt("resolutions do not fit the frequencies and contacts"));
        }
        for (std::size_t f = 0; f < frequencies; ++f) {
            const auto first = resolution.begin() + static_cast<std::ptrdiff_t>(f * contacts);
            p.resolution.emplace_back(first, first + static_cast<std::ptrdiff_t>(contacts));
        }
    }
    NITCAD_IO_TRY(fields, r.records(k::fields));
    for (const Reader& f : fields) {
        NITCAD_IO_TRY(ac, read_ac_fields(f));
        if (ac.contact >= contacts) {
            return std::unexpected(corrupt("small-signal fields of a contact that is not there"));
        }
        p.fields.push_back(std::move(ac));
    }
    return p;
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

namespace {

namespace report_key {
constexpr std::uint16_t counts = 1, min_angle = 2, counts_3d = 3, quality_3d = 4;  // 3D: optional
}  // namespace report_key

std::vector<std::int64_t> wide(std::span<const std::int32_t> v) { return {v.begin(), v.end()}; }

// Each value into [0, top] as an int32, or false.
bool narrow(const std::vector<std::int64_t>& v, std::vector<std::int32_t>& out) {
    for (const std::int64_t x : v) {
        if (x < 0 || x > std::numeric_limits<std::int32_t>::max()) return false;
        out.push_back(static_cast<std::int32_t>(x));
    }
    return true;
}

}  // namespace

std::expected<void, base::Error> save_generated_mesh(const std::filesystem::path& path,
                                                     const mesh::GeneratedMesh& m) {
    Writer w = write_mesh(m.mesh);
    w.i64s(mesh_key::cells, wide(m.simplices.cells));
    w.i64s(mesh_key::facets, wide(m.simplices.facets));
    w.i64s(mesh_key::facet_patch, wide(m.simplices.facet_patch));
    w.i64s(mesh_key::node_region, wide(m.node_region));
    Writer rep;
    const mesh::MeshReport& r = m.report;
    rep.i64s(report_key::counts,
             std::vector<std::int64_t>{static_cast<std::int64_t>(r.nodes), static_cast<std::int64_t>(r.triangles),
                                       static_cast<std::int64_t>(r.interface_nodes),
                                       static_cast<std::int64_t>(r.refinement_points),
                                       static_cast<std::int64_t>(r.segment_splits),
                                       static_cast<std::int64_t>(r.exempt_triangles)});
    rep.f64(report_key::min_angle, r.min_angle_deg);
    if (m.simplices.dimension == 3) {
        rep.i64s(report_key::counts_3d,
                 std::vector<std::int64_t>{static_cast<std::int64_t>(r.tetrahedra), static_cast<std::int64_t>(r.facet_splits),
                                           static_cast<std::int64_t>(r.exempt_tetrahedra), static_cast<std::int64_t>(r.slivers),
                                           static_cast<std::int64_t>(r.flat_dropped),
                                           static_cast<std::int64_t>(r.band_facet_points)});
        rep.f64s(report_key::quality_3d, std::vector<double>{r.max_radius_edge, r.min_dihedral_deg});
    }
    w.record(mesh_key::report, rep);
    if (m.certificate) w.record(mesh_key::certificate, detail::write_certificate(*m.certificate));
    return write_file(path, Content::mesh, w);
}

std::expected<mesh::GeneratedMesh, base::Error> load_generated_mesh(const std::filesystem::path& path) {
    NITCAD_IO_TRY(bytes, read_file(path, Content::mesh));
    NITCAD_IO_TRY(r, Reader::parse(bytes));
    if (!r.has(mesh_key::cells) || !r.has(mesh_key::node_region) || !r.has(mesh_key::report)) {
        return std::unexpected(invalid("'" + path.string() + "' is a plain mesh, not a generated one"));
    }
    NITCAD_IO_TRY(m, read_mesh(r));
    NITCAD_IO_TRY(cells, r.i64s(mesh_key::cells));
    NITCAD_IO_TRY(facets, r.i64s(mesh_key::facets));
    NITCAD_IO_TRY(facet_patch, r.i64s(mesh_key::facet_patch));
    NITCAD_IO_TRY(regions, r.i64s(mesh_key::node_region));
    NITCAD_IO_TRY(rep, r.record(mesh_key::report));
    NITCAD_IO_TRY(counts, rep.i64s(report_key::counts));
    NITCAD_IO_TRY(min_angle, rep.f64(report_key::min_angle));
    if (counts.size() != 6 || regions.size() != m.node_count()) {
        return std::unexpected(corrupt("generated mesh sizes"));
    }
    mesh::SimplexMesh s;
    s.dimension = m.dimension();
    s.points.assign(m.points().begin(), m.points().end());
    std::vector<std::int32_t> node_region;
    if (!narrow(cells, s.cells) || !narrow(facets, s.facets) || !narrow(facet_patch, s.facet_patch) ||
        !narrow(regions, node_region)) {
        return std::unexpected(corrupt("generated mesh index out of range"));
    }
    for (const mesh::BoundaryPatch& b : m.boundary()) s.patch_names.push_back(b.name);
    auto dual = mesh::make_dual_mesh(s);
    if (!dual || mesh::graph_digest(*dual) != mesh::graph_digest(m)) {
        return std::unexpected(corrupt("the simplices do not give the stored mesh"));
    }
    std::optional<mesh::InterfaceCertificate> certificate;
    if (r.has(mesh_key::certificate)) {
        NITCAD_IO_TRY(cr, r.record(mesh_key::certificate));
        NITCAD_IO_TRY(c, detail::read_certificate(cr));
        if (!mesh::certifies(c, m, node_region)) {
            return std::unexpected(corrupt("the certificate does not certify the stored mesh and regions"));
        }
        certificate = c;
    }
    mesh::MeshReport report{static_cast<std::size_t>(counts[0]), static_cast<std::size_t>(counts[1]),
                                  static_cast<std::size_t>(counts[2]), static_cast<std::size_t>(counts[3]),
                                  static_cast<std::size_t>(counts[4]), static_cast<std::size_t>(counts[5]),
                                  min_angle};
    if (rep.has(report_key::counts_3d)) {
        NITCAD_IO_TRY(c3, rep.i64s(report_key::counts_3d));
        NITCAD_IO_TRY(q3, rep.f64s(report_key::quality_3d));
        if (c3.size() != 6 || q3.size() != 2) return std::unexpected(corrupt("generated mesh 3D report sizes"));
        report.tetrahedra = static_cast<std::size_t>(c3[0]);
        report.facet_splits = static_cast<std::size_t>(c3[1]);
        report.exempt_tetrahedra = static_cast<std::size_t>(c3[2]);
        report.slivers = static_cast<std::size_t>(c3[3]);
        report.flat_dropped = static_cast<std::size_t>(c3[4]);
        report.band_facet_points = static_cast<std::size_t>(c3[5]);
        report.max_radius_edge = q3[0];
        report.min_dihedral_deg = q3[1];
    }
    return mesh::GeneratedMesh{std::move(m), std::move(s), std::move(node_region), std::move(certificate), report};
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
    w.record(checkpoint_key::provenance, write_provenance(run.provenance));
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
    if (r.has(checkpoint_key::provenance)) {
        NITCAD_IO_TRY(p, r.record(checkpoint_key::provenance));
        NITCAD_IO_TRY(provenance, read_provenance(p));
        c.provenance = std::move(provenance);
    }
    return c;
}

std::expected<void, base::Error> save_transient(const std::filesystem::path& path,
                                                const results::Transient& t) {
    namespace k = transient_key;
    Writer w;
    w.record(k::run, write_run(t.run));
    for (const results::TimePoint& p : t.points) w.record(k::point, write_time_point(p));
    for (const results::TransientSnapshot& s : t.snapshots) w.record(k::snapshot, write_snapshot(s));
    w.u64(k::rejected, t.rejected_steps);
    w.u64(k::retraced, t.retraced_steps);
    w.u64(k::lagged, t.lagged_steps);
    if (t.stopped) w.record(k::stopped, write_error(*t.stopped));
    if (t.unfinished) w.record(k::unfinished, write_convergence(*t.unfinished));
    return write_file(path, Content::transient, w);
}

std::expected<results::Transient, base::Error> load_transient(const std::filesystem::path& path) {
    namespace k = transient_key;
    NITCAD_IO_TRY(bytes, read_file(path, Content::transient));
    NITCAD_IO_TRY(r, Reader::parse(bytes));
    results::Transient t;
    NITCAD_IO_TRY(run_record, r.record(k::run));
    NITCAD_IO_TRY(run, read_run(run_record));
    t.run = std::move(run);
    NITCAD_IO_TRY(points, r.records(k::point));
    for (const Reader& p : points) {
        NITCAD_IO_TRY(point, read_time_point(p));
        t.points.push_back(std::move(point));
    }
    NITCAD_IO_TRY(snapshots, r.records(k::snapshot));
    for (const Reader& s : snapshots) {
        NITCAD_IO_TRY(snapshot, read_snapshot(s, t.points.size()));
        t.snapshots.push_back(std::move(snapshot));
    }
    NITCAD_IO_TRY(rejected, r.u64(k::rejected));
    NITCAD_IO_TRY(retraced, r.u64(k::retraced));
    NITCAD_IO_TRY(lagged, r.u64(k::lagged));
    NITCAD_IO_TRY(rejected_n, size_of(rejected, "rejected step count"));
    NITCAD_IO_TRY(retraced_n, size_of(retraced, "retraced step count"));
    NITCAD_IO_TRY(lagged_n, size_of(lagged, "lagged step count"));
    t.rejected_steps = rejected_n;
    t.retraced_steps = retraced_n;
    t.lagged_steps = lagged_n;
    if (r.has(k::stopped)) {
        NITCAD_IO_TRY(e, r.record(k::stopped));
        NITCAD_IO_TRY(stopped, read_error(e));
        t.stopped = std::move(stopped);
    }
    if (r.has(k::unfinished)) {
        NITCAD_IO_TRY(c, r.record(k::unfinished));
        NITCAD_IO_TRY(unfinished, read_convergence(c));
        t.unfinished = std::move(unfinished);
    }
    return t;
}

std::expected<void, base::Error> save_small_signal(const std::filesystem::path& path,
                                                   const results::SmallSignal& s) {
    namespace k = small_signal_key;
    Writer w;
    w.record(k::run, write_run(s.run));
    w.u64(k::contacts, s.contacts);
    w.f64s(k::frequencies, s.frequency_Hz);
    for (const results::SmallSignalPoint& p : s.points) w.record(k::point, write_ac_point(p));
    if (s.stopped) w.record(k::stopped, write_error(*s.stopped));
    if (s.unfinished) w.record(k::unfinished, write_convergence(*s.unfinished));
    return write_file(path, Content::small_signal, w);
}

std::expected<results::SmallSignal, base::Error> load_small_signal(
    const std::filesystem::path& path) {
    namespace k = small_signal_key;
    NITCAD_IO_TRY(bytes, read_file(path, Content::small_signal));
    NITCAD_IO_TRY(r, Reader::parse(bytes));
    results::SmallSignal s;
    NITCAD_IO_TRY(run_record, r.record(k::run));
    NITCAD_IO_TRY(run, read_run(run_record));
    s.run = std::move(run);
    NITCAD_IO_TRY(contacts, r.u64(k::contacts));
    NITCAD_IO_TRY(frequencies, r.f64s_or_empty(k::frequencies));
    // A file's arrays are bounded by its size (far below 2^40 values), so these bounds keep
    // frequencies x contacts^2 well inside 64 bits.
    if (contacts > (std::uint64_t{1} << 16) || frequencies.size() > (std::size_t{1} << 30)) {
        return std::unexpected(corrupt("small-signal contact or frequency count"));
    }
    s.contacts = static_cast<std::size_t>(contacts);
    s.frequency_Hz = std::move(frequencies);
    NITCAD_IO_TRY(points, r.records(k::point));
    for (const Reader& p : points) {
        NITCAD_IO_TRY(point, read_ac_point(p, s.frequency_Hz.size(), s.contacts));
        s.points.push_back(std::move(point));
    }
    if (r.has(k::stopped)) {
        NITCAD_IO_TRY(e, r.record(k::stopped));
        NITCAD_IO_TRY(stopped, read_error(e));
        s.stopped = std::move(stopped);
    }
    if (r.has(k::unfinished)) {
        NITCAD_IO_TRY(c, r.record(k::unfinished));
        NITCAD_IO_TRY(unfinished, read_convergence(c));
        s.unfinished = std::move(unfinished);
    }
    return s;
}

}  // namespace NiTCAD::io
