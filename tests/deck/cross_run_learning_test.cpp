// AC19's family gate (Unit 17b, slice 3; DECISIONS.md Unit 17 AC19, ID8-ID10, Unit 17b LB12): 8 tensor
// MOSFETs as decks (AR6's tensor MOSFET with L_g, t_ox and N_A varied), each adapted with the
// intelligent policy, adjoint-weighted, to 1e-2 on I_D at V_g 0.5 and 1 V, its history saved; a model
// trained on 6 by `nitcad learn` and the 2 held out run again with it. On each held-out device the
// learned first mesh needs at least one pass fewer than ID2's a priori field and ends with no more
// DOFs, both certified; the trees beat ridge regression on held-out rank correlation. Manual
// (`[.ac19]`, Release): ten adaptive MOSFET runs, about 19 minutes. Fails as frozen, measured 2026-11-13
// (DECISIONS.md Unit 17b, AC19): m6 learned 4 passes and 12 000 DOFs against ID2's 5 and 10 476 (the DOF
// clause unmet; m7 met, 4 and 10 260 against 5 and 15 561); trees 0.399 against ridge 0.161: a limitation
// (ARCHITECTURE.md 11.3), the criterion unchanged.
#include <catch2/catch_test_macros.hpp>

#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/io/learn_files.hpp"
#include "NiTCAD/results/adapt.hpp"
#include "NiTCAD/solve/adapt_learn.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::must;
namespace fs = std::filesystem;

namespace {

struct Member {
    int lg_nm;
    double tox_nm;
    double na;
};

// The family: 6 trained on, the last 2 held out (inside the training ranges: interpolation).
const std::vector<Member> family{{400, 4.0, 1e17}, {500, 5.0, 2e17}, {600, 5.0, 1e17}, {700, 6.0, 5e16},
                                 {800, 4.0, 2e17}, {500, 3.0, 5e16}, {700, 5.0, 1.5e17}, {600, 4.0, 7e16}};
constexpr std::size_t first_held_out = 6;

std::string mosfet_deck(const Member& m, const std::string& extra) {
    const int L = 600 + m.lg_nm;  // source and drain 300 nm each
    char text[2048];
    std::snprintf(text, sizeof text, R"(
[mesh]
type = tensor
x = linspace(0, %dnm, %d)
y = linspace(0, 200nm, 11)

[region silicon]
material = silicon

[doping]
profile = acceptor constant %g
profile = donor gaussian 1e19 axis=y center=0 width=50nm x=0:300nm lateral=14.142135624nm
profile = donor gaussian 1e19 axis=y center=0 width=50nm x=%dnm:%dnm lateral=14.142135624nm

[contact source]
kind = ohmic
patch = y_min
box = 0 300.5nm  0 0

[contact drain]
kind = ohmic
patch = y_min
box = %gnm %dnm  0 0

[contact body]
kind = ohmic
patch = y_max

[contact gate]
kind = gate
patch = y_min
box = 301nm %dnm  0 0
oxide_thickness = %gnm

[solver]
max_cuts = 20

[sweep]
contact drain = 0.05
contact gate = 0 : 1 : 0.5

[adapt]
goal = current drain
tolerance = 1e-2
at = 1, 2
max_passes = 10
weighting = adjoint
policy = intelligent
)",
                  L, L / 50 + 1, m.na, 300 + m.lg_nm, L, 299.5 + m.lg_nm, L, 299 + m.lg_nm, m.tox_nm);
    return text + extra;
}

struct Outcome {
    std::size_t passes = 0;
    std::uint64_t dofs = 0;
    std::string status, stop, predictor;
};

Outcome run(const fs::path& dir, const std::string& name, const Member& m, const std::string& extra) {
    const deck::Plan plan = must(mosfet_deck(m, extra), dir);
    const auto r = deck::run(plan);
    if (!r) FAIL(name << ": " << r.error().message);
    const results::AdaptReport& a = *r->adapt;
    Outcome o{a.passes.size(), a.passes.back().dofs, a.status, a.stop, a.predictor_status};
    WARN("AC19 " << name << " (L_g " << m.lg_nm << " nm, t_ox " << m.tox_nm << " nm, N_A " << m.na << "): " << o.passes
                 << " passes, " << o.dofs << " DOFs, " << o.status << ", stop " << o.stop
                 << (o.predictor.empty() ? "" : ", predictor " + o.predictor));
    return o;
}

}  // namespace

TEST_CASE("AC19: the learned first mesh on held-out MOSFETs saves a pass at no more DOFs; trees beat ridge",
          "[.ac19][!shouldfail]") {
    const fs::path dir = fs::temp_directory_path() / "nitcad_ac19";
    fs::create_directories(dir);
    std::vector<Outcome> a_priori;
    std::vector<fs::path> histories;
    for (std::size_t k = 0; k < family.size(); ++k) {
        const std::string name = "m" + std::to_string(k);
        a_priori.push_back(run(dir, name + " (a priori)", family[k], "\n[output]\nadapt_history = " + name + ".nithist\n"));
        histories.push_back(dir / (name + ".nithist"));
    }
    auto man = io::make_manifest(histories);
    REQUIRE(man.has_value());
    REQUIRE(io::save_manifest(dir / "family.manifest", *man).has_value());
    auto set = io::build_training_set(*man);
    REQUIRE(set.has_value());
    // The held-out devices by their files' digests (the manifest sorts by digest).
    std::vector<std::size_t> held;
    for (std::size_t k = first_held_out; k < family.size(); ++k) {
        const std::uint64_t d = io::file_digest(histories[k]).value();
        for (std::size_t e = 0; e < man->entries.size(); ++e) {
            if (man->entries[e].digest == d) held.push_back(e);
        }
    }
    REQUIRE(held.size() == 2);
    auto model = solve::train_predictor(*set, held);
    REQUIRE(model.has_value());
    REQUIRE(io::save_learned_model(dir / "family.nitlearn", *model).has_value());
    WARN("AC19 model: " << model->training_samples << " training samples; held out " << model->metrics.samples
                        << " samples, rank correlation trees " << model->metrics.rank_correlation << ", ridge "
                        << model->metrics.ridge_rank_correlation << "; first pass recall " << model->metrics.recall
                        << ", precision " << model->metrics.precision);
    CHECK(model->metrics.rank_correlation > model->metrics.ridge_rank_correlation);

    for (std::size_t k = first_held_out; k < family.size(); ++k) {
        const Outcome learned = run(dir, "m" + std::to_string(k) + " (learned)", family[k], "predictor = family.nitlearn\n");
        const Outcome& prior = a_priori[k];
        CHECK(learned.predictor == "used");
        CHECK(learned.status == "certified");
        CHECK(prior.status == "certified");
        CHECK(learned.passes + 1 <= prior.passes);
        CHECK(learned.dofs <= prior.dofs);
    }
}
