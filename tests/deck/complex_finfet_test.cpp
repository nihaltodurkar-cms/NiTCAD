// FC9 and FC10 (DECISIONS.md Unit 53; owner rulings SQ1, SQ2, 2026-10-30): the two extruded FinFET
// decks solved. Hidden, hours in Release with PARDISO (NITCAD_MKL the path of mkl_rt.3.dll); run by
// tag.
// - [.fc9] FC9: decks/finfet_3d_reduced.deck (generated) against a coarse tensor mesh of
//   finfet_3d.deck's layout and its 2x refined reference. For I_D at every V_G above the current
//   resolution, V_T, max g_m and the swing: |Q_gen - Q_ref| <= 1.5 |Q_tensor - Q_ref| + f, f the
//   reference's Richardson estimate |Q_ref - Q_tensor| / (2^2 - 1) (the box method's second order);
//   and the generated mesh at most the coarse tensor mesh's nodes. FQ2's pre-stated failure path:
//   the generated mesh has about 30 times those nodes (SQ1), so the case is [!shouldfail] with its
//   measurement, the threshold unchanged.
// - [.fc_digests]: the two decks' mesh digests (SQ9).
// - [.fc10] FC10: decks/finfet_3d_shaped.deck meshes, is certified, and its sweep converges at
//   every point; nodes, slivers, mesh and solve times printed (no accuracy gate).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/analysis/dc.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/linalg/mkl_runtime.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;

namespace {

std::string mkl_path() {
    char* value = nullptr;
    std::size_t length = 0;
    const bool given = _dupenv_s(&value, &length, "NITCAD_MKL") == 0 && value != nullptr;
    const std::string mkl = given ? value : "";
    std::free(value);
    return mkl;
}

std::string read_deck(const char* name) {
    std::ifstream in(std::filesystem::path(NITCAD_DECKS_DIR) / name);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

void replace_once(std::string& text, const std::string& from, const std::string& to) {
    const auto at = text.find(from);
    INFO("deck text: " << from);
    REQUIRE(at != std::string::npos);
    REQUIRE(text.find(from, at + 1) == std::string::npos);
    text.replace(at, from.size(), to);
}

std::string pardiso(const std::string& mkl) { return "backend = pardiso\nmkl = " + mkl + "\nthreads = 6\n"; }

// finfet_3d.deck's mesh lines: the coarse mesh halves every linspace's cells (its straddles as the
// deck's; two interfaces 0.5 nm apart allow at most 1/3 nm), the reference is the deck's lines with
// the straddles halved, so the reference is the coarse mesh refined twice everywhere but at the
// straddled pairs, which narrow by half from the deck's.
std::string tensor_text(bool reference, const std::string& mkl) {
    std::string t = read_deck("finfet_3d.deck");
    replace_once(t, "# backend = pardiso\n", pardiso(mkl));
    if (reference) {
        for (const char* axis : {"x = 30nm", "x = 55nm", "y = 20nm"}) {
            replace_once(t, std::string("straddle ") + axis + " 0.5nm", std::string("straddle ") + axis + " 0.25nm");
        }
        for (const char* axis : {"y = 50nm", "y = 50.5nm", "z = 1.5nm", "z = 2nm", "z = 12nm", "z = 12.5nm"}) {
            replace_once(t, std::string("straddle ") + axis + " 0.2nm", std::string("straddle ") + axis + " 0.1nm");
        }
        return t;
    }
    replace_once(t, "x = linspace(0, 20nm, 9) + linspace(20nm, 65nm, 37) + linspace(65nm, 85nm, 9)",
                 "x = linspace(0, 20nm, 5) + linspace(20nm, 65nm, 19) + linspace(65nm, 85nm, 5)");
    replace_once(t, "y = linspace(0, 20nm, 5) + linspace(20nm, 46nm, 14) + linspace(46nm, 50nm, 9) + linspace(50nm, 52nm, 5)",
                 "y = linspace(0, 20nm, 3) + linspace(20nm, 46nm, 7) + linspace(46nm, 50nm, 5) + linspace(50nm, 52nm, 3)");
    replace_once(t,
                 "z = linspace(0, 2nm, 5) + linspace(2nm, 3nm, 3) + linspace(3nm, 11nm, 9) + linspace(11nm, 12nm, 3) + "
                 "linspace(12nm, 14nm, 5)",
                 "z = linspace(0, 2nm, 3) + linspace(2nm, 3nm, 2) + linspace(3nm, 11nm, 5) + linspace(11nm, 12nm, 2) + "
                 "linspace(12nm, 14nm, 3)");
    return t;
}

std::size_t contact_index(const deck::Plan& p, const std::string& name) {
    std::size_t c = 0;
    while (c < p.device.contacts().size() && p.device.contacts()[c].name != name) ++c;
    INFO("contact " << name);
    REQUIRE(c < p.device.contacts().size());
    return c;
}

struct Run {
    std::size_t nodes = 0;
    analysis::Curve transfer;
    double build_s = 0.0, solve_s = 0.0;
    std::size_t converged = 0, requested = 0;
};

Run run(const std::string& label, const std::string& text, const char* gate) {
    const auto t0 = std::chrono::steady_clock::now();
    const deck::Plan p = deck_test::must(text, NITCAD_DECKS_DIR);
    const auto t1 = std::chrono::steady_clock::now();
    auto s = solve::sweep_bias(p.device, p.points, p.options);
    const auto t2 = std::chrono::steady_clock::now();
    if (!s) FAIL(label << ": " << s.error().message);
    Run r;
    r.nodes = p.device.mesh().node_count();
    r.build_s = std::chrono::duration<double>(t1 - t0).count();
    r.solve_s = std::chrono::duration<double>(t2 - t1).count();
    r.converged = s->points.size();
    r.requested = p.points.size();
    std::size_t substeps = 0;
    for (const auto& b : s->points) substeps += b.bias_substeps;
    std::printf("%s: %zu nodes, build %.0f s, sweep %.0f s, %zu of %zu points, %zu bias sub-steps%s%s\n",
                label.c_str(), r.nodes, r.build_s, r.solve_s, r.converged, r.requested, substeps,
                s->stopped ? ", stopped: " : "", s->stopped ? s->stopped->message.c_str() : "");
    if (p.generated) {
        const auto& g = p.generated->report;
        std::printf("%s: %zu tetrahedra, %zu slivers, min dihedral %.2f deg, mesh digest %016llx\n", label.c_str(),
                    g.tetrahedra, g.slivers, g.min_dihedral_deg,
                    static_cast<unsigned long long>(mesh::graph_digest(p.device.mesh())));
    }
    auto c = analysis::current_curve(*s, contact_index(p, gate), contact_index(p, "drain"));
    if (!c) FAIL(label << ": " << c.error().message);
    r.transfer = std::move(*c);
    for (std::size_t k = 0; k < r.transfer.x.size(); ++k) {
        std::printf("%s V_G %.2f V: I_D %.6e A (resolution %.2e)\n", label.c_str(), r.transfer.x[k],
                    r.transfer.y[k], r.transfer.resolution.empty() ? 0.0 : r.transfer.resolution[k]);
    }
    std::fflush(stdout);
    return r;
}

// FC9's rule for one quantity; prints it and returns whether it holds.
bool within(const char* what, double gen, double tensor, double ref) {
    const double d = std::abs(tensor - ref);
    const double bound = 1.5 * d + d / 3.0;
    const bool ok = std::abs(gen - ref) <= bound;
    std::printf("FC9 %s: generated %.6e, coarse tensor %.6e, reference %.6e; |gen - ref| %.3e, bound %.3e%s\n", what,
                gen, tensor, ref, std::abs(gen - ref), bound, ok ? "" : "  FAILS");
    return ok;
}

// FC9's comparison of the deck `name` (meshed by a generator, its gate contact "gate") with the coarse
// tensor mesh and its 2x reference: every point converged, the rule for every quantity, and at most the
// coarse tensor mesh's nodes. UC8 (Unit 54) is the same comparison for the swept deck.
void compare_to_tensor(const char* name) {
    const std::string mkl = mkl_path();
    if (mkl.empty()) SKIP("NITCAD_MKL not set");
    if (auto loaded = linalg::load_mkl_runtime(mkl); !loaded) FAIL(loaded.error().message);
    const Run tensor = run("FC9 coarse tensor", tensor_text(false, mkl), "gate_top");
    const Run ref = run("FC9 tensor reference", tensor_text(true, mkl), "gate_top");
    std::string generated = read_deck(name);
    replace_once(generated, "[solver]\n", "[solver]\n" + pardiso(mkl));
    const Run gen = run(std::string("FC9 generated (") + name + ")", generated, "gate");
    CHECK(tensor.converged == tensor.requested);
    CHECK(ref.converged == ref.requested);
    CHECK(gen.converged == gen.requested);
    REQUIRE(gen.transfer.x == ref.transfer.x);
    REQUIRE(tensor.transfer.x == ref.transfer.x);
    bool all = true;
    for (std::size_t k = 0; k < ref.transfer.x.size(); ++k) {
        const auto above = [k](const Run& r) {
            return r.transfer.resolution.empty() || std::abs(r.transfer.y[k]) > r.transfer.resolution[k];
        };
        if (!above(ref) || !above(tensor) || !above(gen)) continue;
        const std::string what = "I_D at " + std::to_string(ref.transfer.x[k]) + " V";
        all = within(what.c_str(), gen.transfer.y[k], tensor.transfer.y[k], ref.transfer.y[k]) && all;
    }
    const auto extract = [](const Run& r, auto f) {
        const auto e = f(r.transfer);
        if (!e) FAIL(e.error().message);
        return e->value;
    };
    const auto vt = [](const analysis::Curve& c) { return analysis::threshold_max_gm(c, 0.05); };
    const auto gm = [](const analysis::Curve& c) { return analysis::max_transconductance(c); };
    const auto ss = [](const analysis::Curve& c) { return analysis::subthreshold_swing(c); };
    all = within("V_T (max g_m) [V]", extract(gen, vt), extract(tensor, vt), extract(ref, vt)) && all;
    all = within("max g_m [A/V]", extract(gen, gm), extract(tensor, gm), extract(ref, gm)) && all;
    all = within("swing [mV/dec]", extract(gen, ss), extract(tensor, ss), extract(ref, ss)) && all;
    CHECK(all);
    std::printf("FC9 nodes: generated %zu, coarse tensor %zu, reference %zu\n", gen.nodes, tensor.nodes, ref.nodes);
    CHECK(gen.nodes <= tensor.nodes);
}

// The shaped deck `name` solved: every point converges; on-current, V_T and swing printed.
Run shaped_sweep(const char* label, const char* name) {
    const std::string mkl = mkl_path();
    if (mkl.empty()) SKIP("NITCAD_MKL not set");
    if (auto loaded = linalg::load_mkl_runtime(mkl); !loaded) FAIL(loaded.error().message);
    std::string shaped = read_deck(name);
    replace_once(shaped, "[solver]\n", "[solver]\n" + pardiso(mkl));
    const Run r = run(label, shaped, "gate");
    CHECK(r.converged == r.requested);
    const auto vt = analysis::threshold_max_gm(r.transfer, 0.05);
    const auto ss = analysis::subthreshold_swing(r.transfer);
    std::printf("%s: I_D(%.2f V) %.6e A, V_T (max g_m) %.4f V, swing %.2f mV/dec\n", label, r.transfer.x.back(),
                r.transfer.y.back(), vt ? vt->value : std::nan(""), ss ? ss->value : std::nan(""));
    std::fflush(stdout);
    return r;
}

}  // namespace

TEST_CASE("complex FinFET: the reduced deck against a coarse tensor mesh and its 2x reference (FC9)",
          "[.fc9][!shouldfail]") {
    compare_to_tensor("finfet_3d_reduced.deck");
}

// UC8 (Unit 54, the unit's gate): FC9's comparison, unchanged, for the swept deck with structured layers.
// Measured (DECISIONS.md, Unit 54 slice 5): 8,300 nodes, every point converged; V_T, max g_m and the swing within
// their bounds, I_D at 0.8 V +1.0% (the coarse tensor +1.2%); I_D outside its bound at 0, 0.05, 0.1, 0.4 and 0.5 V.
// [!shouldfail] for an owner ruling, the rule unchanged.
TEST_CASE("complex FinFET: the swept deck against a coarse tensor mesh and its 2x reference (UC8)", "[.uc8][!shouldfail]") {
    compare_to_tensor("finfet_3d_swept.deck");
}

TEST_CASE("complex FinFET: the shaped deck's sweep converges at every point (FC10)", "[.fc10]") {
    shaped_sweep("FC10 shaped", "finfet_3d_shaped.deck");
}

// UC9 (Unit 54): the shaped deck as a swept deck; certified (the build refuses otherwise), every point
// converges; nodes, times, on-current, V_T and swing printed beside FC10's run (no accuracy gate).
TEST_CASE("complex FinFET: the shaped swept deck's sweep converges at every point (UC9)", "[.uc9]") {
    shaped_sweep("UC9 shaped swept", "finfet_3d_shaped_swept.deck");
}

// The two decks' meshes pinned (SQ9): the graph digests measured on this branch after SQ6 and equal
// to a build without the FC12 speedups and the review fixes (DECISIONS.md Unit 53, slice 6). Hidden:
// about 7 minutes each in Release.
TEST_CASE("complex FinFET: the two decks' meshes are pinned", "[.fc_digests]") {
    for (const auto& [name, digest] : {std::pair{"finfet_3d_reduced.deck", 0xd5a99737270489bfULL},
                                       std::pair{"finfet_3d_shaped.deck", 0x719a79d9189eba8aULL}}) {
        const deck::Plan p = deck_test::must(read_deck(name), NITCAD_DECKS_DIR);
        INFO(name);
        CHECK(mesh::graph_digest(p.device.mesh()) == digest);
    }
}
