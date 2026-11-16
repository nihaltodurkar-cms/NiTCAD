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

#include "NiTCAD/assemble/fields.hpp"
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
            REQUIRE(dofs.fields()[k].unknown);
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
