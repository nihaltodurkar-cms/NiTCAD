// New device (Unit 26n slice 1, DT4): the library's families, a deck's main dimensions, a thumbnail of its regions,
// and the copy a Create writes.
#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "NiTCAD/app/new_device.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;
namespace devices = NiTCAD::app::devices;

namespace {

const fs::path devices_dir = fs::path(NITCAD_EXAMPLES_DIR) / "devices";

std::string read_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

}  // namespace

TEST_CASE("new device: the library's families, in file-name order, with their dimension", "[app][new_device]") {
    const auto f = devices::list_families(devices_dir);
    REQUIRE(f.size() == 29);
    std::vector<std::string> names;
    std::vector<int> dims;
    for (const auto& x : f) {
        names.push_back(x.example.path.filename().string());
        dims.push_back(x.dimension);
    }
    CHECK(names.front() == "algan_gan_hemt_2d.deck");
    CHECK(f.front().example.title == "AlGaN/GaN HEMT (2D)");
    CHECK(dims == std::vector<int>{2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 1, 2, 3, 1, 1, 2, 3, 1, 2, 2, 3, 2, 2, 3, 1, 2, 2, 2});
    // A deck that does not build is listed, dimension 0.
    const fs::path dir = app_test::scratch("new_device_families");
    app_test::write(dir, "broken.deck", "# Broken\n# It does not build.\n[mesh]\ntype = tensor\n");
    const auto b = devices::list_families(dir);
    REQUIRE(b.size() == 1);
    CHECK(b[0].dimension == 0);
    CHECK(b[0].example.title == "Broken");
}

TEST_CASE("new device: a deck's main dimensions are its first group's rows", "[app][new_device]") {
    CHECK(devices::main_dimensions(read_file(devices_dir / "finfet_bulk_3d.deck")) ==
          std::vector<std::string>{"Lg", "Wfin", "Hfin", "tox", "Nsd"});
    CHECK(devices::main_dimensions("[mesh]\ntype = tensor\n").empty());
    CHECK(devices::main_dimensions("[parameters]\na = 1nm\nb = 2nm\n# group: Other\nc = 3nm\n") ==
          std::vector<std::string>{"a", "b"});
}

TEST_CASE("new device: a thumbnail of the regions, cut at the middle of z", "[app][new_device]") {
    const auto m = devices::thumbnail(read_file(devices_dir / "mosfet_2d.deck"));
    CHECK(m.error.empty());
    CHECK(m.x0 == 0.0);
    CHECK(m.y0 == 0.0);
    CHECK(std::abs(m.x1 - 2e-4) < 1e-12);
    CHECK(std::abs(m.y1 - 1e-4) < 1e-12);
    REQUIRE(m.boxes.size() == 2);
    CHECK(m.boxes[0].material == "silicon");  // the background: the whole domain
    CHECK(m.boxes[1].material == "sio2");
    CHECK(std::abs(m.boxes[1].y0 - 0.99e-4) < 1e-12);
    CHECK(std::abs(m.boxes[1].y1 - 1e-4) < 1e-12);
    // 3D: the fin's box crosses the middle of z; every box cut there.
    const auto f = devices::thumbnail(read_file(devices_dir / "finfet_bulk_3d.deck"));
    CHECK(f.error.empty());
    CHECK(f.boxes.size() == 6);  // air, insulator, gate block, gate oxide, substrate, fin
    // 1D: one bar of unit height.
    const auto p = devices::thumbnail(read_file(devices_dir / "pn_diode_1d.deck"));
    REQUIRE(p.boxes.size() == 1);
    CHECK(p.y1 - p.y0 == 1.0);
    // A deck that does not build: its message, no boxes.
    const auto e = devices::thumbnail("[mesh]\ntype = tensor\n");
    CHECK_FALSE(e.error.empty());
    CHECK(e.boxes.empty());
}

TEST_CASE("new device: Create writes the deck under the typed name, never over a file", "[app][new_device]") {
    const fs::path dir = app_test::scratch("new_device_create") / "not" / "yet";
    const std::string text = "# T\n[mesh]\ntype = tensor\n";
    const auto a = devices::create_device(text, dir, "my_fet");
    REQUIRE(a);
    CHECK(*a == dir / "my_fet.deck");
    CHECK(read_file(*a) == text);
    const auto b = devices::create_device(text, dir, "other.deck");
    REQUIRE(b);
    CHECK(*b == dir / "other.deck");
    // An existing file stays as it is.
    std::ofstream(*a, std::ios::binary | std::ios::app) << "# mine\n";
    const auto again = devices::create_device(text, dir, "my_fet");
    REQUIRE_FALSE(again);
    CHECK(again.error() == "my_fet.deck is already in " + dir.string() + ": choose another name");
    CHECK(read_file(*a) == text + "# mine\n");
    for (const char* bad : {"a/b", R"(a\b)", "x:y", "", "  ", "q?", "a*"}) {
        CAPTURE(bad);
        CHECK_FALSE(devices::create_device(text, dir, bad));
    }
}

TEST_CASE("new device: Create names the run's outputs after the new deck, so two devices of one family never share them", "[app][new_device]") {
    const fs::path dir = app_test::scratch("new_device_outputs");
    const std::string text = "# T\n[parameters]\nL = 1um\n\n[output]\nsweep = mosfet_2d.nitsweep\ndevice = mosfet_2d.nitdev  # the device\ntable = mosfet_2d.csv\n";
    const auto a = devices::create_device(text, dir, "fet_a");
    REQUIRE(a);
    CHECK(read_file(*a) == "# T\n[parameters]\nL = 1um\n\n[output]\nsweep = fet_a.nitsweep\ndevice = fet_a.nitdev  # the device\ntable = fet_a.csv\n");
}
