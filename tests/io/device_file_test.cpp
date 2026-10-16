// Device files (Unit 31, io/device_file.hpp): every field of a description round-trips bit for
// bit; a loaded description is validated by Device::create; damaged and foreign files are refused.
#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "NiTCAD/io/device_file.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/mesh/interfaces.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "container.hpp"
#include "record.hpp"

using namespace NiTCAD;
using base::ErrorCode;
using device::ContactKind;
using device::Device;
using device::DeviceDescription;
namespace fs = std::filesystem;

namespace {

fs::path scratch(const std::string& name) {
    const fs::path dir = fs::temp_directory_path() / "nitcad_io_test" / name;
    fs::remove_all(dir);
    fs::create_directories(dir);
    return dir;
}

std::vector<char> bytes_of(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

std::vector<mesh::NodeId> patch(const mesh::Mesh& m, const char* name) {
    const mesh::BoundaryPatch* b = m.find_boundary(name);
    REQUIRE(b != nullptr);
    return {b->nodes.begin(), b->nodes.end()};
}

// A 2D device using every part of a description: silicon, an AlGaAs region behind a thermionic
// interface, an oxide on top with an electrode, fixed charge, trap levels and bands and surface
// recombination at the oxide, graded doping, tensor cells and two thermal contacts.
Device rich_device() {
    const std::vector<double> x{0, 0.2e-4, 0.4e-4, 0.55e-4, 0.65e-4, 0.8e-4, 1e-4};
    const std::vector<double> y{0, 0.1e-4, 0.2e-4, 0.29e-4, 0.31e-4, 0.4e-4};
    auto m = *mesh::make_tensor_grid(x, y);
    auto cells = *mesh::TensorCells::create(m, x, y);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, 0);
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        if (p[1] > 0.3e-4) {
            region[i] = 2;
            continue;
        }
        if (p[0] > 0.6e-4) region[i] = 1;
        donors[i] = 1e15 * (1.0 + p[0] * 1e4);
        acceptors[i] = p[0] < 0.5e-4 ? 3.1e17 : 0.0;
    }
    std::vector<mesh::NodeId> top, left, right;
    for (const mesh::NodeId k : patch(m, "y_max")) top.push_back(k);
    for (const mesh::NodeId k : patch(m, "x_min")) {
        if (region[static_cast<std::size_t>(k)] != 2) left.push_back(k);
    }
    for (const mesh::NodeId k : patch(m, "x_max")) {
        if (region[static_cast<std::size_t>(k)] != 2) right.push_back(k);
    }
    auto oxide_params = physics::silicon_dioxide_parameters;
    oxide_params.eps_r = 3.9000000000000004;  // one ulp off: must survive
    device::Interface oxide{.region_a = "silicon", .region_b = "oxide",
                            .fixed_charge_cm2 = -2.5e10,
                            .recombination_velocity_n_cm_s = 1e4,
                            .recombination_velocity_p_cm_s = 2e4};
    oxide.traps.levels = {{physics::TrapType::donor, 1e10, -0.1, 2e-15, 3e-16},
                          {physics::TrapType::acceptor, 2e10, 0.2, 1e-15, 1e-15}};
    oxide.traps.bands = {{physics::TrapType::acceptor, 1e11, 0.0, 0.4, 1e-16, 2e-16}};
    oxide.traps.thermal_velocity_n_cm_s = 2.3e7;
    auto d = Device::create(
        {.mesh = std::move(m),
         .temperature_K = 301.5,
         .regions = {{"silicon", physics::silicon()},
                     {"algaas", *physics::Semiconductor::create(*physics::algaas_parameters(0.3))},
                     {"oxide", *physics::Insulator::create(oxide_params)}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", ContactKind::ohmic, std::move(left)},
                      {"cathode", ContactKind::ohmic, std::move(right)},
                      {"gate", ContactKind::electrode, std::move(top), {},
                       {.kind = device::GateElectrode::metal, .work_function_eV = 4.61}}},
         .interfaces = {oxide,
                        {"silicon", "algaas", device::InterfaceTransport::thermionic_emission}},
         .cells = std::move(cells),
         .thermal_contacts = {{"sink", "x_min", device::ThermalContactKind::isothermal, 310.0, 0.0},
                              {"package", "x_max", device::ThermalContactKind::resistance, 290.0,
                               1e-4}}});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

// A 1D silicon diode with a lumped gate on a graph mesh, whose two regions need a P4 certificate.
Device certified_device() {
    std::vector<mesh::Point> p;
    std::vector<mesh::Edge> e;
    for (int i = 0; i < 5; ++i) p.push_back({1e-4 * i, 0, 0});
    for (int i = 0; i < 4; ++i) e.push_back({i, i + 1, 1e-4, 1.0});
    auto m = *mesh::Mesh::from_parts(1, p, std::vector<double>(5, 1e-4), e,
                                     {{"x_min", {0}, {1.0}}, {"x_max", {4}, {1.0}}},
                                     mesh::MeshOrigin::graph);
    const std::vector<device::RegionId> region{0, 0, 1, 1, 1};
    const std::vector<mesh::InterfacePoint1D> at{{1.5e-4, 0, 1}};
    auto certificate = *mesh::certify_interfaces_1d(m, region, at);
    auto d = Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"left", physics::silicon()}, {"right", physics::silicon()}},
         .node_region = region,
         .donors = {1e16, 1e16, 0, 0, 0},
         .acceptors = {0, 0, 1e17, 1e17, 1e17},
         .contacts = {{"a", ContactKind::ohmic, {0}},
                      {"g", ContactKind::gate, {4},
                       {.boundary = "x_max", .oxide_thickness_cm = 1e-6,
                        .electrode = device::GateElectrode::p_poly, .fixed_charge_cm2 = 3e10}}},
         .interface_certificate = certificate});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

bool same_bits(std::span<const double> a, std::span<const double> b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::bit_cast<std::uint64_t>(a[i]) != std::bit_cast<std::uint64_t>(b[i])) return false;
    }
    return true;
}

void check_same(const Device& a, const Device& b) {
    CHECK(io::mesh_digest(a.mesh()) == io::mesh_digest(b.mesh()));
    CHECK(a.mesh().origin() == b.mesh().origin());
    CHECK(std::bit_cast<std::uint64_t>(a.temperature_K()) ==
          std::bit_cast<std::uint64_t>(b.temperature_K()));
    REQUIRE(a.regions().size() == b.regions().size());
    for (std::size_t r = 0; r < a.regions().size(); ++r) {
        CHECK(a.regions()[r].name == b.regions()[r].name);
        CHECK(a.regions()[r].material.index() == b.regions()[r].material.index());
        if (const auto* s = std::get_if<physics::Semiconductor>(&a.regions()[r].material)) {
            CHECK(s->parameters() ==
                  std::get<physics::Semiconductor>(b.regions()[r].material).parameters());
        } else {
            CHECK(std::get<physics::Insulator>(a.regions()[r].material).parameters() ==
                  std::get<physics::Insulator>(b.regions()[r].material).parameters());
        }
    }
    CHECK(std::vector(a.node_region().begin(), a.node_region().end()) ==
          std::vector(b.node_region().begin(), b.node_region().end()));
    CHECK(same_bits(a.donors(), b.donors()));
    CHECK(same_bits(a.acceptors(), b.acceptors()));
    REQUIRE(a.contacts().size() == b.contacts().size());
    for (std::size_t c = 0; c < a.contacts().size(); ++c) {
        const auto& x = a.contacts()[c];
        const auto& y = b.contacts()[c];
        CHECK(x.name == y.name);
        CHECK(x.kind == y.kind);
        CHECK(x.nodes == y.nodes);
        CHECK(x.gate.boundary == y.gate.boundary);
        CHECK(x.gate.oxide_thickness_cm == y.gate.oxide_thickness_cm);
        CHECK(x.gate.oxide_relative_permittivity == y.gate.oxide_relative_permittivity);
        CHECK(x.gate.electrode == y.gate.electrode);
        CHECK(x.gate.work_function_eV == y.gate.work_function_eV);
        CHECK(x.gate.fixed_charge_cm2 == y.gate.fixed_charge_cm2);
        CHECK(x.electrode.kind == y.electrode.kind);
        CHECK(x.electrode.work_function_eV == y.electrode.work_function_eV);
    }
    REQUIRE(a.interfaces().size() == b.interfaces().size());
    for (std::size_t f = 0; f < a.interfaces().size(); ++f) {
        const auto& x = a.interfaces()[f];
        const auto& y = b.interfaces()[f];
        CHECK(x.region_a == y.region_a);
        CHECK(x.region_b == y.region_b);
        CHECK(x.transport == y.transport);
        CHECK(x.fixed_charge_cm2 == y.fixed_charge_cm2);
        CHECK(x.traps == y.traps);
        CHECK(x.recombination_velocity_n_cm_s == y.recombination_velocity_n_cm_s);
        CHECK(x.recombination_velocity_p_cm_s == y.recombination_velocity_p_cm_s);
    }
    REQUIRE((a.cells() == nullptr) == (b.cells() == nullptr));
    if (a.cells()) {
        for (int k = 0; k < a.cells()->dimension(); ++k) {
            CHECK(same_bits(a.cells()->axis(k), b.cells()->axis(k)));
        }
    }
    REQUIRE(a.thermal_contacts().size() == b.thermal_contacts().size());
    for (std::size_t t = 0; t < a.thermal_contacts().size(); ++t) {
        const auto& x = a.thermal_contacts()[t];
        const auto& y = b.thermal_contacts()[t];
        CHECK(x.name == y.name);
        CHECK(x.boundary == y.boundary);
        CHECK(x.kind == y.kind);
        CHECK(x.temperature_K == y.temperature_K);
        CHECK(x.resistance_K_cm2_W == y.resistance_K_cm2_W);
    }
    REQUIRE((a.interface_certificate() == nullptr) == (b.interface_certificate() == nullptr));
    if (a.interface_certificate()) {
        const auto& x = *a.interface_certificate();
        const auto& y = *b.interface_certificate();
        CHECK(x.mesh_digest == y.mesh_digest);
        CHECK(x.region_digest == y.region_digest);
        CHECK(x.interfaces_digest == y.interfaces_digest);
        CHECK(x.report.cross_faces == y.report.cross_faces);
        CHECK(x.report.corner_faces == y.report.corner_faces);
        CHECK(x.report.max_flat_deviation == y.report.max_flat_deviation);
        CHECK(x.report.max_corner_deviation == y.report.max_corner_deviation);
    }
}

// Writes a device payload by hand (the container of files.hpp around `w`).
void write_raw(const fs::path& p, const io::detail::Writer& w) {
    REQUIRE(io::detail::write_file(p, io::Content::device, w).has_value());
}

// The minimal valid 1D device payload, with the silicon region's parameter record built by `params`.
template <class F>
io::detail::Writer minimal_payload(F&& params, std::uint64_t region_kind = 0) {
    using io::detail::Writer;
    const auto m = *mesh::make_tensor_grid(std::vector<double>{0, 1e-4, 2e-4});
    Writer w;
    w.record(1, io::detail::write_mesh(m));
    w.f64(2, 300.0);
    Writer region;
    region.text(1, "si");
    region.u64(2, region_kind);
    Writer p;
    params(p);
    region.record(3, p);
    w.record(3, region);
    w.i64s(4, std::vector<std::int64_t>{0, 0, 0});
    w.f64s(5, std::vector<double>{1e16, 1e16, 1e16});
    w.f64s(6, std::vector<double>{0, 0, 0});
    Writer c;
    c.text(1, "a");
    c.u64(2, 0);
    c.i64s(3, std::vector<std::int64_t>{0});
    w.record(7, c);
    return w;
}

// Legacy silicon's parameters by the device file's keys (device_file.cpp), as a file written
// before Unit 35 has them, skipping `skip`.
void silicon_keys(io::detail::Writer& p, int skip = 0) {
    const auto& s = physics::legacy_silicon_parameters;
    const std::vector<std::pair<int, double>> v{
        {1, s.eps_r}, {2, s.Eg0_eV}, {3, s.varshni_alpha_eV_per_K}, {4, s.varshni_beta_K},
        {5, s.Nc300}, {6, s.Nv300}, {7, s.electron_affinity_eV},
        {10, s.electron_mobility.mu_min}, {11, s.electron_mobility.mu_max},
        {12, s.electron_mobility.N_ref}, {13, s.electron_mobility.alpha},
        {14, s.electron_mobility.T_exponent}, {15, s.hole_mobility.mu_min},
        {16, s.hole_mobility.mu_max}, {17, s.hole_mobility.N_ref}, {18, s.hole_mobility.alpha},
        {19, s.hole_mobility.T_exponent}, {20, s.lifetime.tau_n0}, {21, s.lifetime.tau_p0},
        {22, s.lifetime.N_ref}, {23, s.auger.Cn}, {24, s.auger.Cp},
        {25, s.bandgap_narrowing.E0_eV}, {26, s.bandgap_narrowing.N0},
        {27, s.electron_saturation.v_sat_cm_s}, {28, s.electron_saturation.beta},
        {29, s.hole_saturation.v_sat_cm_s}, {30, s.hole_saturation.beta},
        {32, s.ionization.donor_eV}, {33, s.ionization.acceptor_eV},
        {34, s.ionization.donor_degeneracy}, {35, s.ionization.acceptor_degeneracy}};
    for (const auto& [k, x] : v) {
        if (k != skip) p.f64(static_cast<std::uint16_t>(k), x);
    }
}

}  // namespace

TEST_CASE("device file: every field round-trips bit for bit", "[io][device_file]") {
    const fs::path dir = scratch("device_round_trip");
    for (const Device& d : {rich_device(), certified_device()}) {
        const fs::path first = dir / "first.nitdev";
        const fs::path second = dir / "second.nitdev";
        REQUIRE(io::save_device(first, d).has_value());
        const auto loaded = io::load_device(first);
        if (!loaded) FAIL(loaded.error().message);
        check_same(d, *loaded);
        // Saving the loaded device reproduces the file byte for byte: nothing is lost or reordered.
        REQUIRE(io::save_device(second, *loaded).has_value());
        CHECK(bytes_of(first) == bytes_of(second));
        CHECK_FALSE(fs::exists(dir / "first.nitdev.tmp"));
    }
}

TEST_CASE("device file: foreign and damaged files are refused", "[io][device_file]") {
    const fs::path dir = scratch("device_refused");
    const Device d = certified_device();
    const fs::path file = dir / "d.nitdev";
    REQUIRE(io::save_device(file, d).has_value());

    SECTION("a mesh file is not a device file") {
        const fs::path mfile = dir / "m.nitmesh";
        REQUIRE(io::save_mesh(mfile, d.mesh()).has_value());
        const auto r = io::load_device(mfile);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().code == ErrorCode::invalid_input);
        CHECK(r.error().message.find("content kind") != std::string::npos);
        CHECK_FALSE(io::load_mesh(file).has_value());
    }
    SECTION("a flipped payload byte") {
        auto b = bytes_of(file);
        b[40] = static_cast<char>(b[40] ^ 0x10);
        std::ofstream(file, std::ios::binary | std::ios::trunc).write(b.data(), static_cast<std::streamsize>(b.size()));
        const auto r = io::load_device(file);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().message.find("checksum") != std::string::npos);
    }
    SECTION("a missing file") {
        const auto r = io::load_device(dir / "absent.nitdev");
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().code == ErrorCode::invalid_input);
    }
}

TEST_CASE("device file: a loaded description is validated and versioned by key", "[io][device_file]") {
    const fs::path dir = scratch("device_keys");
    const fs::path file = dir / "d.nitdev";

    SECTION("the hand-built minimal file (no Unit 35 keys) loads as legacy silicon") {
        write_raw(file, minimal_payload([](io::detail::Writer& p) { silicon_keys(p); }));
        const auto r = io::load_device(file);
        if (!r) FAIL(r.error().message);
        const auto& s = std::get<physics::Semiconductor>(r->regions()[0].material).parameters();
        auto expected = physics::legacy_silicon_parameters;
        // Keys the file leaves out take the parameter struct's defaults (IO2): no impact
        // ionization, tunnelling or thermal data, no radiative or Richardson terms.
        expected.impact_ionization = {};
        expected.band_to_band = {};
        expected.thermal = {};
        CHECK(s == expected);
    }
    SECTION("an unknown key is skipped") {
        write_raw(file, minimal_payload([](io::detail::Writer& p) {
                      silicon_keys(p);
                      p.f64(999, 1.0);
                  }));
        CHECK(io::load_device(file).has_value());
    }
    SECTION("an absent required parameter is refused by the material") {
        write_raw(file, minimal_payload([](io::detail::Writer& p) { silicon_keys(p, 1); }));
        const auto r = io::load_device(file);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().code == ErrorCode::invalid_input);
        CHECK(r.error().message.find("eps_r") != std::string::npos);
    }
    SECTION("an invalid parameter is refused by the material") {
        write_raw(file, minimal_payload([](io::detail::Writer& p) {
                      silicon_keys(p, 2);
                      p.f64(2, -1.0);
                  }));
        REQUIRE_FALSE(io::load_device(file).has_value());
    }
    SECTION("an unknown region kind") {
        write_raw(file, minimal_payload([](io::detail::Writer& p) { silicon_keys(p); }, 7));
        const auto r = io::load_device(file);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().message.find("region kind") != std::string::npos);
    }
    SECTION("a description Device::create rejects") {
        auto bad = minimal_payload([](io::detail::Writer& p) { silicon_keys(p); });
        io::detail::Writer c;
        c.text(1, "a");
        c.u64(2, 0);
        c.i64s(3, std::vector<std::int64_t>{1});  // a second contact named "a", on an inner node
        bad.record(7, c);
        write_raw(file, bad);
        const auto r = io::load_device(file);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().message.find("repeated") != std::string::npos);
    }
}
