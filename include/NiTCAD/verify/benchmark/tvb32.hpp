// TVB-32, the consistency suite C1-C12 on real solves of the three TVR-1 decks (DECISIONS Unit 55 11.6, TVR1-Q5,
// the TVB-32 rulings and RC-Q5 of 11.13). The decks are embedded at build time (verify_decks.hpp, generated from
// examples/pn_diode_1d.deck, decks/schottky_diode_1d.deck, examples/moscap_cv.deck) and read by deck::parse, outputs
// cleared. Every check whose inputs the run supplies is applied at every operating point: C1-C5 and C12 at each DC
// point, C1 and C3 on the equilibrium, C7 and C8 at each C-V point and at the zero-bias AC point added through the
// API. C6, C9, C10 and C11 need inputs no deck produces (listed in the benchmark's notes).
#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/results/small_signal.hpp"
#include "NiTCAD/results/solution.hpp"
#include "NiTCAD/results/transient.hpp"
#include "NiTCAD/verify/consistency/checks.hpp"

namespace NiTCAD::verify::benchmark {

struct Tvb32Entry {
    std::string deck;  // "pn_diode_1d", "schottky_diode_1d", "moscap_cv"
    std::string at;    // "equilibrium", "bias (0.1) V", "added AC point, bias (0) V"
    bool zero_bias = false;  // every bias exactly 0 and one thermal temperature (the C1 predicates)
    std::size_t contacts = 0;
    consistency::CheckResult result;
};

struct Tvb32Run {
    std::vector<Tvb32Entry> entries;
    std::optional<std::string> run_error;
};

[[nodiscard]] Tvb32Run tvb32_entries();

// The TVR-4 devices (DECISIONS Unit 55 16.44), TVB-32's full tier: each benchmark's own device at the coarsest level
// of its ladder, at its own operating points plus zero bias. C1-C5 and C12 at each DC point, C1 and C3 on the
// equilibrium (devices with a semiconductor region), C7 and C8 at each small-signal point, C6 on each transient, C9
// on TVB-18's diode.
struct Tvb32Collector {
    std::vector<Tvb32Entry>& out;
    std::string deck;  // the benchmark and its device, e.g. "TVB-17 diode N = 64"
    std::size_t contacts = 0;
    void add(const std::string& at, bool zero, consistency::CheckResult r);
    void dc(const results::BiasPoint& p, const mesh::Mesh& m, const consistency::DeviceData& d,
            const std::string& at);
    void ac(const results::SmallSignal& run, const mesh::Mesh& m, const consistency::DeviceData& d,
            const std::string& prefix);
    void transient(const results::Transient& run, const mesh::Mesh& m, const consistency::DeviceData& d,
                   const std::string& at);
    // C1 and C3 on solve_equilibrium with these models; false (error set) when it does not solve.
    bool equilibrium(const device::Device& device, const assemble::PhysicsModels& models, std::string& error);
};
[[nodiscard]] bool at_zero_bias(const results::BiasPoint& p);
[[nodiscard]] std::string bias_label(const results::BiasPoint& p);

[[nodiscard]] Tvb32Run tvb05_consistency();
[[nodiscard]] Tvb32Run tvb06_consistency();
[[nodiscard]] Tvb32Run tvb10_consistency();
[[nodiscard]] Tvb32Run tvb13_consistency();
[[nodiscard]] Tvb32Run tvb15_consistency();
[[nodiscard]] Tvb32Run tvb16_consistency();
[[nodiscard]] Tvb32Run tvb17_consistency();
[[nodiscard]] Tvb32Run tvb18_consistency();
[[nodiscard]] Tvb32Run tvb19_consistency();
[[nodiscard]] Tvb32Run tvb20_consistency();  // with step 9 (16.46)
[[nodiscard]] Tvb32Run tvb27_consistency();
[[nodiscard]] Tvb32Run tvb28_consistency();
// All of them, in that order; a run error stops the list at that benchmark.
[[nodiscard]] Tvb32Run tvb32_tvr4_entries();

}  // namespace NiTCAD::verify::benchmark
