// Fields and the unknown map (Unit 40 slice 1; DECISIONS.md FE1, FG7): the drift-diffusion
// system's map reproduces the node-interleaved layout (psi, n, p, then the temperature rise with
// the electrothermal model), and every unknown's field and node can be asked for.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/assemble/band_edges.hpp"
#include "NiTCAD/assemble/charges.hpp"
#include "NiTCAD/assemble/currents.hpp"
#include "NiTCAD/assemble/derived.hpp"
#include "NiTCAD/assemble/fields.hpp"
#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/assemble/rates.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/assemble/system.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD;
using assemble::FieldId;

namespace {

device::Device resistor(bool thermal) {
    std::vector<double> x(11);
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = 1e-5 * static_cast<double>(i);
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    std::vector<device::ThermalContact> sinks;
    if (thermal) {
        sinks = {{"left", "x_min", device::ThermalContactKind::isothermal, 300.0, 0.0},
                 {"right", "x_max", device::ThermalContactKind::isothermal, 300.0, 0.0}};
    }
    return *device::Device::create({.mesh = std::move(m),
                                    .temperature_K = 300.0,
                                    .regions = {{"silicon", physics::silicon()}},
                                    .node_region = std::vector<device::RegionId>(n, 0),
                                    .donors = std::vector<double>(n, 1e16),
                                    .acceptors = std::vector<double>(n, 0.0),
                                    .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                                                 {"right", device::ContactKind::ohmic, std::move(right)}},
                                    .thermal_contacts = std::move(sinks)});
}

}  // namespace

TEST_CASE("fields: the drift-diffusion map is the node-interleaved layout") {
    for (const bool thermal : {false, true}) {
        CAPTURE(thermal);
        const device::Device d = resistor(thermal);
        const auto scaling = *assemble::make_scaling(d);
        const auto system =
            *assemble::System::create(d, scaling, {.electrothermal = thermal});
        const assemble::DofMap& dofs = system.dofs();
        const std::size_t m = thermal ? 4 : 3;
        REQUIRE(dofs.block() == m);
        REQUIRE(system.stride() == m);
        REQUIRE(dofs.node_count() == system.node_count());
        REQUIRE(dofs.unknowns() == system.unknowns());
        REQUIRE(dofs.has(FieldId::temperature) == thermal);
        const FieldId order[] = {FieldId::potential, FieldId::electrons, FieldId::holes,
                                 FieldId::temperature};
        for (std::size_t k = 0; k < m; ++k) {
            REQUIRE(dofs.slot(order[k]) == k);
            REQUIRE(dofs.fields()[k].id == order[k]);
            REQUIRE(dofs.fields()[k].kind == assemble::FieldKind::unknown);
            REQUIRE(dofs.fields()[k].location == assemble::FieldLocation::node);
        }
        for (std::size_t i = 0; i < system.node_count(); ++i) {
            for (std::size_t k = 0; k < m; ++k) {
                const std::size_t u = dofs.index(order[k], i);
                REQUIRE(u == m * i + k);
                REQUIRE(dofs.field_of(u).id == order[k]);
                REQUIRE(dofs.node_of(u) == i);
            }
        }
    }
}

TEST_CASE("fields: the descriptors") {
    REQUIRE(assemble::potential_field.domain == assemble::FieldDomain::all);
    REQUIRE(assemble::electron_field.domain == assemble::FieldDomain::semiconductor);
    REQUIRE(assemble::hole_field.domain == assemble::FieldDomain::semiconductor);
    REQUIRE(assemble::temperature_field.domain == assemble::FieldDomain::all);
    const assemble::Field two[] = {assemble::hole_field, assemble::potential_field};
    const auto map = assemble::DofMap::node_interleaved(two, 5);
    REQUIRE(map.unknowns() == 10);
    REQUIRE(map.slot(FieldId::holes) == 0);
    REQUIRE(map.slot(FieldId::potential) == 1);
    REQUIRE(!map.slot(FieldId::electrons).has_value());
    REQUIRE(map.index(FieldId::potential, 3) == 7);
    REQUIRE(map.field_of(7).id == FieldId::potential);
    REQUIRE(map.node_of(7) == 3);
}

TEST_CASE("derived fields: the table (Unit 57 GE4)") {
    using assemble::DiagnosticId;
    const auto table = assemble::diagnostics();
    REQUIRE(table.size() == 6);
    std::vector<FieldId> seen;
    for (std::size_t r = 0; r < table.size(); ++r) {
        const assemble::Diagnostic& d = table[r];
        CAPTURE(d.name);
        REQUIRE(static_cast<std::size_t>(d.id) == r);
        REQUIRE(&assemble::diagnostic(d.id) == &d);
        REQUIRE(!d.fields.empty());
        REQUIRE(d.evaluate != nullptr);
        const auto location = d.id == DiagnosticId::edge_currents ||
                                      d.id == DiagnosticId::edge_flux_terms
                                  ? assemble::FieldLocation::edge
                              : d.id == DiagnosticId::trap_occupancies
                                  ? assemble::FieldLocation::interface_edge
                                  : assemble::FieldLocation::node;
        for (std::size_t k = 0; k < d.fields.size(); ++k) {
            const assemble::Field& f = d.fields[k];
            CAPTURE(f.name);
            REQUIRE(f.kind == assemble::FieldKind::derived);
            REQUIRE(f.location == location);
            REQUIRE(!f.units.empty());
            REQUIRE(assemble::position(d, f.id) == k);
            REQUIRE(std::ranges::find(seen, f.id) == seen.end());
            seen.push_back(f.id);
        }
    }
    REQUIRE(seen.size() == 21);
    // No unknown is a derived field.
    for (const FieldId u : {FieldId::potential, FieldId::electrons, FieldId::holes,
                            FieldId::temperature, FieldId::electron_temperature,
                            FieldId::hole_temperature}) {
        REQUIRE(std::ranges::find(seen, u) == seen.end());
    }
}

TEST_CASE("derived fields: the values are the topic functions' (Unit 57 GE4)") {
    using assemble::DiagnosticId;
    struct Case {
        const char* name;
        bool thermal;
        assemble::PhysicsModels models;
    };
    const Case cases[] = {
        {"drift-diffusion", false, {}},
        {"electrothermal", true, {.electrothermal = true}},
        {"energy balance + impact", false, {.impact_ionization = true, .energy_balance = true}},
    };
    for (const Case& c : cases) {
        CAPTURE(c.name);
        const device::Device d = resistor(c.thermal);
        const auto scaling = *assemble::make_scaling(d);
        const auto system = *assemble::System::create(d, scaling, c.models);
        const std::size_t nodes = system.node_count(), m = system.stride();
        // A state off equilibrium: psi rising across the bar, the electrons raised in the
        // interior, and every extra field (a rise) varying from node to node.
        std::vector<double> psi(nodes);
        for (std::size_t i = 0; i < nodes; ++i) psi[i] = 4.0 + 0.3 * static_cast<double>(i);
        std::vector<double> x = system.state_from_potential(psi);
        for (std::size_t i = 1; i + 1 < nodes; ++i) {
            x[m * i + 1] *= 1.0 + 0.05 * static_cast<double>(i);
            for (std::size_t k = 3; k < m; ++k) x[m * i + k] = 0.01 * static_cast<double>(i + k);
        }

        auto eval = [&](DiagnosticId id) {
            return assemble::diagnostic(id).evaluate(system, x);
        };
        const auto& bands_row = assemble::diagnostic(DiagnosticId::bands);
        auto bands = eval(DiagnosticId::bands);
        const assemble::BandEdges b = band_edges(system, x);
        REQUIRE(take(bands, bands_row, FieldId::conduction_band) == b.conduction);
        REQUIRE(take(bands, bands_row, FieldId::valence_band) == b.valence);
        REQUIRE(take(bands, bands_row, FieldId::electron_fermi) == b.electron_fermi);
        REQUIRE(take(bands, bands_row, FieldId::hole_fermi) == b.hole_fermi);

        const auto& rates_row = assemble::diagnostic(DiagnosticId::rates);
        auto r = eval(DiagnosticId::rates);
        REQUIRE(r.size() == rates_row.fields.size());
        std::vector<double> net_n(nodes, 0.0), net_p(nodes, 0.0);
        std::size_t present = 0;
        for (const assemble::MechanismRate& mr : rates(system, x)) {
            for (std::size_t i = 0; i < nodes; ++i) {
                net_n[i] += mr.electron[i];
                net_p[i] += mr.hole[i];
            }
            using M = assemble::RateMechanism;
            const auto at = [&](FieldId f) { return r[assemble::position(rates_row, f)]; };
            switch (mr.mechanism) {
                case M::srh: REQUIRE(at(FieldId::srh_rate) == mr.electron); break;
                case M::auger: REQUIRE(at(FieldId::auger_rate) == mr.electron); break;
                case M::radiative: REQUIRE(at(FieldId::radiative_rate) == mr.electron); break;
                case M::impact: REQUIRE(at(FieldId::impact_rate) == mr.electron); break;
                case M::interfaces:
                    REQUIRE(at(FieldId::interface_rate_n) == mr.electron);
                    REQUIRE(at(FieldId::interface_rate_p) == mr.hole);
                    break;
                case M::band_to_band:
                    REQUIRE(at(FieldId::band_to_band_rate_n) == mr.electron);
                    REQUIRE(at(FieldId::band_to_band_rate_p) == mr.hole);
                    break;
            }
            ++present;
        }
        REQUIRE(present >= 3);  // SRH, Auger, radiative (and impact ionization)
        REQUIRE(r[assemble::position(rates_row, FieldId::impact_rate)].empty() ==
                !c.models.impact_ionization);
        REQUIRE(r[assemble::position(rates_row, FieldId::interface_rate_n)].empty());
        REQUIRE(take(r, rates_row, FieldId::net_rate_n) == net_n);
        REQUIRE(take(r, rates_row, FieldId::net_rate_p) == net_p);

        const auto& currents_row = assemble::diagnostic(DiagnosticId::edge_currents);
        auto j = eval(DiagnosticId::edge_currents);
        const auto pairs = edge_currents(system, x);
        const std::vector<double> jn = take(j, currents_row, FieldId::electron_edge_current);
        const std::vector<double> jp = take(j, currents_row, FieldId::hole_edge_current);
        REQUIRE(jn.size() == pairs.size());
        REQUIRE(jp.size() == pairs.size());
        for (std::size_t e = 0; e < pairs.size(); ++e) {
            REQUIRE(jn[e] == pairs[e].first);
            REQUIRE(jp[e] == pairs[e].second);
        }
        REQUIRE(pairs.front().first != pairs.front().second);  // so a swap is seen

        auto flux = eval(DiagnosticId::edge_flux_terms);
        REQUIRE(flux.size() == 1);
        REQUIRE(flux[0] == edge_flux_terms(system, x));
        auto traps = eval(DiagnosticId::trap_occupancies);
        REQUIRE(traps.size() == 1);
        REQUIRE(traps[0] == trap_occupancies(system, x));

        const auto& t_row = assemble::diagnostic(DiagnosticId::temperatures);
        auto t = eval(DiagnosticId::temperatures);
        const double T0 = scaling.temperature_K;
        const auto expected = [&](FieldId rise) {
            std::vector<double> v;
            if (!system.dofs().has(rise)) return v;
            const std::size_t k = *system.dofs().slot(rise);
            for (std::size_t i = 0; i < nodes; ++i) v.push_back(T0 + T0 * x[m * i + k]);
            return v;
        };
        REQUIRE(take(t, t_row, FieldId::lattice_temperature) == expected(FieldId::temperature));
        REQUIRE(take(t, t_row, FieldId::electron_temperature_K) ==
                expected(FieldId::electron_temperature));
        REQUIRE(take(t, t_row, FieldId::hole_temperature_K) ==
                expected(FieldId::hole_temperature));
        REQUIRE(system.dofs().has(FieldId::temperature) == c.models.electrothermal);
        REQUIRE(system.dofs().has(FieldId::electron_temperature) == c.models.energy_balance);
    }
}

namespace {

std::string read_file(const std::filesystem::path& p) {
    std::ifstream in(p, std::ios::binary);
    std::ostringstream s;
    s << in.rdbuf();
    return s.str();
}

// The text with its // and /* */ comments blanked (string literals are not looked into: the
// header has none).
std::string without_comments(const std::string& text) {
    std::string out = text;
    for (std::size_t i = 0; i + 1 < out.size(); ++i) {
        if (out[i] == '/' && out[i + 1] == '/') {
            while (i < out.size() && out[i] != '\n') out[i++] = ' ';
        } else if (out[i] == '/' && out[i + 1] == '*') {
            const std::size_t end = out.find("*/", i + 2);
            const std::size_t last = end == std::string::npos ? out.size() : end + 2;
            for (; i < last; ++i) out[i] = out[i] == '\n' ? '\n' : ' ';
        }
    }
    return out;
}

// The names declared in a class body: every identifier directly followed (spaces aside) by one
// of ( ; = { , ) [ -- the member functions, data members, nested types and their members, and
// the parameters; a type name (followed by <, &, *, :: or a name) is not one.
std::vector<std::string> declared_names(const std::string& body) {
    std::vector<std::string> names;
    const auto ident = [](char c) {
        return std::isalnum(static_cast<unsigned char>(c)) != 0 || c == '_';
    };
    for (std::size_t i = 0; i < body.size();) {
        if (!ident(body[i]) || std::isdigit(static_cast<unsigned char>(body[i]))) {
            ++i;
            continue;
        }
        std::size_t j = i;
        while (j < body.size() && ident(body[j])) ++j;
        std::size_t k = j;
        while (k < body.size() && std::isspace(static_cast<unsigned char>(body[k]))) ++k;
        constexpr std::string_view after = "(;={,)[";
        if (k < body.size() && after.find(body[k]) != std::string_view::npos) {
            names.push_back(body.substr(i, j - i));
        }
        i = j;
    }
    return names;
}

}  // namespace

TEST_CASE("system: the public header names no model and the old type is retired (Unit 57 GG7)") {
    const std::filesystem::path root = NITCAD_SOURCE_DIR;
    REQUIRE_FALSE(std::filesystem::exists(root / "include/NiTCAD/assemble/drift_diffusion.hpp"));
    for (const char* internal : {"src/assemble/terms.hpp", "src/assemble/system_data.hpp",
                                 "src/assemble/composition.hpp"}) {
        CAPTURE(internal);
        CHECK(read_file(root / internal).find("TermAccess") == std::string::npos);
    }
    const std::string header =
        without_comments(read_file(root / "include/NiTCAD/assemble/system.hpp"));
    REQUIRE(header.find("DriftDiffusion") == std::string::npos);
    // The class body, by brace depth.
    const std::size_t start = header.find("class System {");
    REQUIRE(start != std::string::npos);
    std::size_t end = header.find('{', start);
    for (int depth = 0; end < header.size(); ++end) {
        if (header[end] == '{') ++depth;
        if (header[end] == '}' && --depth == 0) break;
    }
    REQUIRE(end < header.size());
    const std::vector<std::string> names = declared_names(header.substr(start, end - start));
    // The parse found the members (a parse that finds nothing would pass the scan).
    for (const char* member : {"create", "dofs", "unknowns", "set_bias", "make_jacobian",
                               "evaluate", "residual", "storage", "small_signal_matrix",
                               "bias_derivative", "bias_rows", "update_size", "apply_update",
                               "stamp_contacts", "state_from_potential", "TimeStep", "data_"}) {
        CAPTURE(member);
        CHECK(std::ranges::find(names, std::string(member)) != names.end());
    }
    for (std::string name : names) {
        CAPTURE(name);
        std::ranges::transform(name, name.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
        for (const char* model : {"thermal", "heat", "tunnel", "path", "energy", "temperature",
                                  "schottky", "gate"}) {
            CHECK(name.find(model) == std::string::npos);
        }
    }
}
