// Shared by the image tests (Unit 25a, VG4) and the band surface's (Unit 26o, BSC8): the WARP renderer, the
// reference comparison, and the small solved MOSFETs. Moved unchanged from image_test.cpp.
//
// The references are tests/render/references/*.png, rendered on WARP. A render matches its reference when no
// more than 0.5% of its pixels differ by more than 24 in any channel. To write new references, run with
// NITCAD_UPDATE_REFERENCES=1, look at every image, and commit them with the reason.
#pragma once

#include <catch2/catch_test_macros.hpp>

#include <cstdlib>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/render/compare.hpp"
#include "NiTCAD/render/image.hpp"
#include "NiTCAD/solve/bias.hpp"

namespace render_test {

using namespace NiTCAD;
using namespace NiTCAD::render;
namespace fs = std::filesystem;

inline const fs::path references = NITCAD_RENDER_REFERENCES;
inline const fs::path output = fs::temp_directory_path() / "nitcad_render_test";

inline Renderer& warp() {
    static Renderer r = [] {
        auto created = Renderer::create({.adapter = Adapter::warp,
#ifndef NDEBUG
                                         .debug_layer = true
#endif
        });
        if (!created) FAIL(created.error().message);
        return std::move(*created);
    }();
    return r;
}

// Draws, saves the render next to the test's output, and compares with the reference (or writes it).
inline void check_image(const Image& image, const std::string& name) {
    fs::create_directories(output);
    if (auto saved = save_png(image, output / (name + ".png")); !saved) FAIL(saved.error().message);
    const fs::path ref = references / (name + ".png");
    char* update = nullptr;
    std::size_t length = 0;
    const bool updating = _dupenv_s(&update, &length, "NITCAD_UPDATE_REFERENCES") == 0 && update &&
                          std::string(update) == "1";
    std::free(update);
    if (updating) {
        fs::create_directories(references);
        REQUIRE(save_png(image, ref).has_value());
        WARN("reference written: " << ref.string());
        return;
    }
    INFO("reference " << ref.string() << ", render " << (output / (name + ".png")).string());
    const auto expected = load_png(ref);
    REQUIRE(expected.has_value());
    const auto d = compare(image, *expected, 24);
    REQUIRE(d.has_value());
    CAPTURE(d->max_channel, d->fraction_over);
    CHECK(d->fraction_over <= 0.005);
}

inline std::vector<double> uniform(double a, double b, int n) {
    std::vector<double> v;
    for (int k = 0; k < n; ++k) v.push_back(k == n - 1 ? b : a + (b - a) * k / (n - 1));
    return v;
}

// A small n-channel MOSFET in 2D or 3D (cm): silicon under an oxide 0.1 um thick on top
// (y in [0.9, 1] um), n+ source and drain at the top corners, a p channel (1e17) between them over a
// lightly doped substrate (1e15); an electrode on the oxide over the channel, ohmic source and drain
// on the pockets' x faces, the body on y_min. Solved at Vg = 1 V, Vd = 0.5 V.
struct Mosfet {
    device::Device device;
    results::BiasPoint point;
    std::vector<results::BiasPoint> sweep;  // every point, for the filmstrip (Unit 25c)
};

inline const Mosfet& mosfet(bool three_d) {
    static auto make = [](bool three) {
        const auto x = uniform(0.0, 2e-4, 21);
        const auto y = uniform(0.0, 1e-4, 21);
        const auto z = three ? uniform(0.0, 0.8e-4, 7) : std::vector<double>{};
        mesh::Mesh m = three ? *mesh::make_tensor_grid(x, y, z) : *mesh::make_tensor_grid(x, y);
        const std::size_t n = m.node_count();
        std::vector<device::RegionId> region(n);
        std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const auto& p = m.points()[i];
            region[i] = p[1] > 0.92e-4 ? 0 : 1;
            if (region[i] == 0) continue;
            const bool top = p[1] > 0.55e-4;
            acceptors[i] = top ? 1e17 : 1e15;
            if (top && (p[0] < 0.55e-4 || p[0] > 1.45e-4)) donors[i] = 2e19;
        }
        std::vector<mesh::NodeId> gate, source, drain;
        for (const mesh::NodeId v : m.find_boundary("y_max")->nodes) {
            const double px = m.points()[static_cast<std::size_t>(v)][0];
            if (px > 0.6e-4 && px < 1.4e-4) gate.push_back(v);
        }
        for (const char* side : {"x_min", "x_max"}) {
            for (const mesh::NodeId v : m.find_boundary(side)->nodes) {
                const double py = m.points()[static_cast<std::size_t>(v)][1];
                if (py > 0.6e-4 && py < 0.92e-4) (side[2] == 'm' && side[3] == 'i' ? source : drain).push_back(v);
            }
        }
        auto body = m.find_boundary("y_min")->nodes;
        device::Contact g{"gate", device::ContactKind::electrode, std::move(gate)};
        g.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.1};
        auto cells = three ? *mesh::TensorCells::create(m, x, y, z) : *mesh::TensorCells::create(m, x, y);
        auto d = device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
             .node_region = std::move(region),
             .donors = std::move(donors),
             .acceptors = std::move(acceptors),
             .contacts = {std::move(g),
                          {"source", device::ContactKind::ohmic, std::move(source)},
                          {"drain", device::ContactKind::ohmic, std::move(drain)},
                          {"body", device::ContactKind::ohmic, std::move(body)}},
             .cells = std::move(cells)});
        if (!d) FAIL(d.error().message);
        const std::vector<std::vector<double>> points{{0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0},
                                                      {1.0, 0.0, 0.0, 0.0}, {1.0, 0.0, 0.25, 0.0},
                                                      {1.0, 0.0, 0.5, 0.0}};
        auto sweep = solve::sweep_bias(*d, points);
        if (!sweep) FAIL(sweep.error().message);
        if (sweep->stopped) FAIL(sweep->stopped->message);
        return Mosfet{std::move(*d), sweep->points.back(), sweep->points};
    };
    static const Mosfet two = make(false);
    static const Mosfet three = make(true);
    return three_d ? three : two;
}

}  // namespace render_test
