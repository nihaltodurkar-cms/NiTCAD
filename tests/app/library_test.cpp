// The device library (Unit 26n slice 1, DT3, DC2-DC4): parametric decks in examples/devices, each titled and
// described, every geometric number and doping level a [parameters] row, its first group "Main"; each builds
// within its node limit, moves by one panel edit per row, and runs at its defaults.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <format>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "NiTCAD/app/dimensions.hpp"
#include "NiTCAD/app/examples.hpp"
#include "NiTCAD/app/run_controller.hpp"
#include "NiTCAD/app/workspace.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;
namespace dims = NiTCAD::app::dims;

namespace {

const fs::path devices_dir = fs::path(NITCAD_EXAMPLES_DIR) / "devices";

std::vector<app::Example> library() { return app::list_examples(devices_dir); }

std::string read_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

std::vector<std::string> lines_of(const std::string& s) {
    std::vector<std::string> out;
    std::size_t at = 0;
    while (true) {
        const auto nl = s.find('\n', at);
        out.push_back(s.substr(at, nl == std::string::npos ? std::string::npos : nl - at));
        if (nl == std::string::npos) return out;
        at = nl + 1;
    }
}

int changed_lines(const std::string& a, const std::string& b) {
    const auto x = lines_of(a), y = lines_of(b);
    if (x.size() != y.size()) return -1;
    int n = 0;
    for (std::size_t k = 0; k < x.size(); ++k) n += x[k] != y[k];
    return n;
}

struct Built {
    bool ok = false;
    std::string error;
    std::size_t nodes = 0;
    int dimension = 0;
};

Built build(const std::string& text) {
    const auto d = deck::parse(text);
    if (!d) return {false, d.error().message};
    const auto plan = deck::build(*d);
    if (!plan) return {false, plan.error().message};
    return {true, {}, plan->device.mesh().node_count(), plan->device.mesh().dimension()};
}

// Keys whose values are not geometry or doping (materials, kinds, work functions, mesh resolution controls).
bool exempt_key(const std::string& k) {
    static const std::vector<std::string> keys{"type", "kind", "material", "electrode", "work_function", "eps_r", "patch",
                                               "richardson_n", "richardson_p", "fixed_charge", "max_size",
                                               "interface_resolution", "layered"};
    return std::ranges::find(keys, k) != keys.end() || k.starts_with("straddle");
}

// The numeric literals of a value that are not 0 and not a linspace count (the last argument before ')').
std::vector<std::string> literal_numbers(const std::string& v) {
    std::vector<std::string> out;
    for (std::size_t i = 0; i < v.size();) {
        const char c = v[i];
        const bool starts = std::isdigit(static_cast<unsigned char>(c)) || (c == '.' && i + 1 < v.size() && std::isdigit(static_cast<unsigned char>(v[i + 1])));
        const char prev = i ? v[i - 1] : ' ';
        if (!starts || std::isalnum(static_cast<unsigned char>(prev)) || prev == '_' || prev == '$' || prev == '.') {
            ++i;
            continue;
        }
        std::size_t j = i;
        while (j < v.size() && (std::isdigit(static_cast<unsigned char>(v[j])) || v[j] == '.')) ++j;
        if (j < v.size() && (v[j] == 'e' || v[j] == 'E') && j + 1 < v.size() &&
            (std::isdigit(static_cast<unsigned char>(v[j + 1])) || v[j + 1] == '-' || v[j + 1] == '+')) {
            j += 2;
            while (j < v.size() && std::isdigit(static_cast<unsigned char>(v[j]))) ++j;
        }
        const std::string num = v.substr(i, j - i);
        while (j < v.size() && std::isalpha(static_cast<unsigned char>(v[j]))) ++j;  // a unit
        std::size_t after = j;
        while (after < v.size() && v[after] == ' ') ++after;
        const bool count = after < v.size() && v[after] == ')';
        if (std::stod(num) != 0.0 && !count) out.push_back(v.substr(i, j - i));
        i = j;
    }
    return out;
}

}  // namespace

TEST_CASE("library: every deck is titled, described, grouped and builds within its node limit", "[app][library]") {
    const auto decks = library();
    REQUIRE_FALSE(decks.empty());
    for (const app::Example& e : decks) {
        CAPTURE(e.path.filename().string());
        CHECK_FALSE(e.title.empty());
        CHECK_FALSE(e.description.empty());
        const std::string text = read_file(e.path);
        const auto rows = dims::read_dimensions(text);
        REQUIRE(rows);
        REQUIRE_FALSE(rows->empty());
        CHECK(rows->front().group == "Main");
        const Built b = build(text);
        INFO(b.error);
        REQUIRE(b.ok);
        std::printf("library: %s %dD, %zu nodes\n", e.path.filename().string().c_str(), b.dimension, b.nodes);
        CHECK(b.nodes < (b.dimension == 3 ? 20000u : 5000u));
    }
}

TEST_CASE("library: every geometric number and doping level is a parameter (DT3)", "[app][library]") {
    for (const app::Example& e : library()) {
        CAPTURE(e.path.filename().string());
        const auto d = deck::parse(read_file(e.path));
        REQUIRE(d);
        for (const deck::Section& s : d->sections) {
            if (s.kind != "mesh" && s.kind != "region" && s.kind != "contact" && s.kind != "doping") continue;
            for (const deck::Entry& en : s.entries) {
                if (exempt_key(en.key)) continue;
                CAPTURE(en.line, en.key, en.value);
                CHECK(literal_numbers(en.value).empty());
            }
        }
    }
}

TEST_CASE("library: one panel edit changes one line, and every length moved by δ still builds with the same node count (DC2, DC3)", "[app][library]") {
    for (const app::Example& e : library()) {
        CAPTURE(e.path.filename().string());
        const std::string text = read_file(e.path);
        const Built base = build(text);
        REQUIRE(base.ok);
        const auto rows = dims::read_dimensions(text);
        REQUIRE(rows);
        for (const dims::Dimension& r : *rows) {
            if (!r.length || r.expr.find('$') != std::string::npos) continue;
            CAPTURE(r.name, r.expr);
            const double step = std::max(1e-8, std::round(0.01 * std::abs(r.value) / 1e-8) * 1e-8);
            const auto out = dims::set_dimension(text, r.name, dims::length_literal(r, r.value + step));
            REQUIRE(out);
            CHECK(changed_lines(text, *out) == 1);
            const Built moved = build(*out);
            INFO(moved.error);
            REQUIRE(moved.ok);
            CHECK(moved.nodes == base.nodes);
        }
    }
}

TEST_CASE("library: the families (DT3)", "[app][library]") {
    std::vector<std::string> names;
    for (const app::Example& e : library()) names.push_back(e.path.filename().string());
    CHECK(names == std::vector<std::string>{
                       "algan_gan_hemt_2d.deck", "avalanche_diode_2d.deck", "back_to_back_diodes_2d.deck", "dual_gate_mosfet_2d.deck",
                       "finfet_bulk_3d.deck", "finfet_soi_3d.deck", "gaa_nanosheet_fet_3d.deck",
                       "gaa_nanowire_fet_3d.deck", "heterojunction_diode_2d.deck", "jfet_2d.deck", "mis_diode_2d.deck",
                       "moscap_1d.deck", "mosfet_2d.deck", "mosfet_3d.deck", "pin_diode_1d.deck", "pn_diode_1d.deck",
                       "pn_diode_2d.deck", "pn_diode_3d.deck", "schottky_diode_1d.deck", "schottky_diode_2d.deck", "sige_hbt_2d.deck",
                       "soi_fd_mosfet_3d.deck", "soi_mosfet_2d.deck", "soi_pd_mosfet_2d.deck", "tri_gate_mosfet_3d.deck",
                       "tunnel_diode_1d.deck", "tunnel_fet_2d.deck", "varactor_diode_2d.deck", "zener_diode_2d.deck"});
}

namespace {

// A library deck run to completion at its defaults; a transfer curve's threshold and swing finite, the swing
// at least 60 mV/dec less 1 % (DC4).
void runs_at_defaults(const app::Example& e, const fs::path& dir) {
    CAPTURE(e.path.filename().string());
    auto copy = app::copy_example(e, dir);
    REQUIRE(copy);
    app::RunController rc;
    const auto t0 = std::chrono::steady_clock::now();
    REQUIRE(rc.start_run(read_file(*copy), *copy));
    rc.wait();
    rc.process();
    INFO((rc.error() ? rc.error()->message : std::string()));
    REQUIRE(rc.state() == app::RunState::done);
    CHECK_FALSE(rc.stopped_early());
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    auto taken = rc.take_results();
    REQUIRE(taken);
    auto w = app::workspace_from_run(std::move(taken->plan), std::move(taken->report));
    REQUIRE(w);
    std::string line = e.path.filename().string() + std::format(": {:.1f} s", seconds);
    double vt = NAN, ss = NAN;
    for (const app::Figure& f : app::extracted_figures(*w)) {
        if (f.name == "Threshold (max gm)") vt = f.value;
        if (f.name == "Subthreshold swing") ss = f.value;
    }
    const std::string stem = e.path.stem().string();
    const bool transistor = stem.find("mosfet") != std::string::npos || stem.find("finfet") != std::string::npos;  // DC4
    if (transistor) {
        line += std::format(", V_th {:.3f} V, SS {:.1f} mV/dec", vt, ss);
        CHECK(std::isfinite(vt));
        CHECK(std::isfinite(ss));
        CHECK(ss >= 60.0 * 0.99);
    }
    std::printf("DC4: %s\n", line.c_str());
}

}  // namespace

TEST_CASE("library: the 1D and 2D decks run to completion at their defaults (DC4)", "[.library_runs]") {
    const fs::path dir = app_test::scratch("library_run");
    for (const app::Example& e : library()) {
        if (build(read_file(e.path)).dimension < 3) runs_at_defaults(e, dir);
    }
}

TEST_CASE("library: the 3D decks run to completion at their defaults (DC4)", "[.library_runs3d]") {
    const fs::path dir = app_test::scratch("library_run_3d");
    for (const app::Example& e : library()) {
        if (build(read_file(e.path)).dimension == 3) runs_at_defaults(e, dir);
    }
}

TEST_CASE("library: the 2D MOSFET at Lg + 1 nm equals the deck written out by hand (DC3)", "[app][library]") {
    const std::string text = read_file(devices_dir / "mosfet_2d.deck");
    const auto rows = dims::read_dimensions(text);
    REQUIRE(rows);
    const auto lg = std::ranges::find_if(*rows, [](const dims::Dimension& r) { return r.name == "Lg"; });
    REQUIRE(lg != rows->end());
    REQUIRE(lg->expr == "1um");
    const auto moved = dims::set_dimension(text, "Lg", "1.001um");
    REQUIRE(moved);
    const std::string by_hand = R"([mesh]
type = tensor
x = linspace(0, 0.5um, 11) + linspace(0.5um, 1.501um, 51) + linspace(1.501um, 2.001um, 11)
y = linspace(0, 0.8um, 17) + linspace(0.8um, 0.99um, 39) + linspace(0.99um, 1um, 6)
straddle y = 0.99um 0.5nm

[region silicon]
material = silicon

[region oxide]
material = sio2
box = 0 2.001um  0.99um 1um

[doping]
profile = acceptor constant 1e17
profile = donor gaussian 1e20 axis=y center=0.99um width=0.1um x=0:0.55um lateral=20nm
profile = donor gaussian 1e20 axis=y center=0.99um width=0.1um x=1.451um:2.001um lateral=20nm

[contact source]
kind = ohmic
patch = x_min
box = 0 0  0.8um 0.985um

[contact drain]
kind = ohmic
patch = x_max
box = 2.001um 2.001um  0.8um 0.985um

[contact gate]
kind = electrode
patch = y_max
box = 0.5um 1.501um  1um 1um
electrode = metal
work_function = 4.1

[contact body]
kind = ohmic
patch = y_min

[equilibrium]
)";
    const auto da = deck::parse(*moved), db = deck::parse(by_hand);
    REQUIRE(da);
    REQUIRE(db);
    const auto pa = deck::build(*da), pb = deck::build(*db);
    REQUIRE(pa);
    INFO((pb ? std::string() : pb.error().message));
    REQUIRE(pb);
    const auto& a = pa->device;
    const auto& b = pb->device;
    const auto close = [](double x, double y) { return std::abs(x - y) <= 1e-12 * std::max({std::abs(x), std::abs(y), 1e-6}); };
    REQUIRE(a.mesh().node_count() == b.mesh().node_count());
    int coords = 0, regions = 0, doping = 0;
    for (std::size_t k = 0; k < a.mesh().node_count(); ++k) {
        for (int c = 0; c < 3; ++c) coords += !close(a.mesh().points()[k][c], b.mesh().points()[k][c]);
        regions += a.node_region()[k] != b.node_region()[k];
        doping += !close(a.donors()[k], b.donors()[k]) || !close(a.acceptors()[k], b.acceptors()[k]);
    }
    CHECK(coords == 0);
    CHECK(regions == 0);
    CHECK(doping == 0);
    for (const auto& c : b.contacts()) {
        CAPTURE(c.name);
        const auto* ca = a.find_contact(c.name);
        REQUIRE(ca);
        CHECK(ca->nodes == c.nodes);
    }
}
